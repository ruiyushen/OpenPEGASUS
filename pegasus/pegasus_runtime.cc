#include "pegasus/pegasus_runtime.h"

#include <seal/seal.h>
#include <seal/util/polyarithsmallmod.h>

namespace gemini {
PegasusRunTime::PegasusRunTime(Parms parms, size_t num_threads)
    : parms_(parms), num_threads_(std::max<size_t>(1, num_threads)) {
  size_t nlevels = parms.nlevels;

  if (parms.enable_repacking) {
    nlevels += (BetterSine::Depth());
  }

  std::string JSON = level2Params(nlevels);

  runtime_ = gemini::RunTime::Create(JSON);
  runtime_->ShowContext(std::cout);
  setUpRuntime();
  setUpSecretKeys();
  setUpPublicKeys();
  setUpFunctors();

  printf("n_{lwe} = %d, n_{lut} = %d, n_{ckks} = %d\n", parms_.lvl0_lattice_dim,
         parms_.lvl1_lattice_dim, parms_.lvl2_lattice_dim);
  printf("KS_base = 2^%d, sk.hamming = %d\n", KS_DC_BASE, SECRET_KEY_HW);
  printf("|msg| < %f, scale = 2^%f, extra_scale = 2^%f, nslots = %d",
         MsgRange(), std::log2(parms_.scale), std::log2(ExtraScaling()),
         parms_.nslots);
  printf(", #thread = %zd\n", num_threads_);
}

PegasusRunTime::PegasusRunTime(Parms parms, size_t num_threads, const std::string &runtime_file) 
 : parms_(parms), num_threads_(std::max<size_t>(1, num_threads)){
    size_t nlevels = parms.nlevels;
  std::string JSON = level2Params(nlevels);
    if (parms.enable_repacking) {
      nlevels += (BetterSine::Depth());
    }

  runtime_ = gemini::RunTime::Create(JSON);
    LoadFromFile(runtime_file, JSON);
    if (!runtime_) {
        throw std::runtime_error("Can't load RunTime");
    }
    runtime_->ShowContext(std::cout);

    setUpRuntime();
    setUpSecretKeys();
    setUpPublicKeys();
    setUpFunctors();

    printf("n_{lwe} = %d, n_{lut} = %d, n_{ckks} = %d\n", parms_.lvl0_lattice_dim,
           parms_.lvl1_lattice_dim, parms_.lvl2_lattice_dim);
    printf("KS_base = 2^%d, sk.hamming = %d\n", KS_DC_BASE, SECRET_KEY_HW);
    printf("|msg| < %f, scale = 2^%f, extra_scale = 2^%f, nslots = %d",
           MsgRange(), std::log2(parms_.scale), std::log2(ExtraScaling()),
           parms_.nslots);
    printf(", #thread = %zd\n", num_threads_);
}

Status PegasusRunTime::SlotsToCoeffs(Ctx *ct, int nct) const {
  ThreadPool pool(num_threads_);
  const size_t work_load = (nct + num_threads_ - 1) / num_threads_;
  for (int w = 0; w < num_threads_; ++w) {
    size_t start = w * work_load;
    size_t end = std::min<size_t>(start + work_load, nct);
    if (end > start) {
      pool.enqueue(
          [&](size_t s, size_t e) {
            for (size_t i = s; i < e; ++i) {
              CHECK_AND_ABORT(SlotsToCoeffs(ct[i]));
            }
          },
          start, end);
    }
  }
  return Status::Ok();
}

Status PegasusRunTime::SlotsToCoeffs(Ctx &out) const {
  size_t level = GetNModuli(out) - 1;
  size_t depth = linearTransformer_->depth();
  if (level < depth) {
    return Status::NotReady("SlotsToCoeffs require more levels");
  }
  // We keep only one moduli after S2C
  runtime_->DropModuli(&out, level - linearTransformer_->s2c_lvl_start());
  Ctx s2c;
  CHECK_STATUS(linearTransformer_->SlotsToCoeffs(out, &s2c));
  out = s2c;
  return Status::Ok();
}

void PegasusRunTime::load_cipher(Ctx &ct, std::ifstream &is) const {
  ct.load(runtime_->SEALRunTime(),is);
}

Status PegasusRunTime::RotateLeft(Ctx &out, size_t offset) const {
  CHECK_STATUS(runtime_->RotateLeft(&out, offset));
  return Status::Ok();
}

Status PegasusRunTime::Add(Ctx &a, Ctx const &b) const {
  CHECK_STATUS(runtime_->Add(&a, b));
  return Status::Ok();
}

Status PegasusRunTime::Sub(Ctx &a, Ctx const &b) const {
  if (&a == &b) {
    return Status::NotReady("Sub itself is not supported");
  }
  CHECK_STATUS(runtime_->Sub(&a, b));
  return Status::Ok();
}

Status PegasusRunTime::Mul(Ctx &a, const Ctx &b) const {
  if (a.parms_id() != b.parms_id()) {
      return Status::ArgumentError("Mul: mismatched parameters in ciphertexts");
  }
  CHECK_STATUS(runtime_->Mul(&a, b));
  return Status::Ok();
}


Status PegasusRunTime::Square(Ctx &a) const {
  CHECK_STATUS(runtime_->Mul(&a, a));
  return Status::Ok();
}

Status PegasusRunTime::RelinThenRescale(Ctx &a) const {
  CHECK_STATUS(runtime_->Relin(&a));
  F64 scale_up = parms_.scale * runtime_->GetModulusPrime(GetNModuli(a) - 1);
  scale_up = std::round(scale_up / a.scale());
  if (scale_up >= 1.) {
    CHECK_STATUS(runtime_->MulScalar(&a, 1., scale_up));
  }
  CHECK_STATUS(runtime_->RescaleNext(&a));
  return Status::Ok();
}

Status PegasusRunTime::ExtraAllCoefficients(const Ctx &in,
                                            std::vector<lwe::Ctx_st> &lwe_ct) {
  std::vector<rlwe::RLWE2LWECt_st> lwe_N_ct(parms_.nslots);
  if (!in.is_ntt_form()) {
    rlwe::SampleExtract(lwe_N_ct.data(), in, extract_indices_,
                        runtime_->SEALRunTime());
  } else {
    auto copy{in};
    rlwe::SwitchNTTForm(copy, runtime_->SEALRunTime());
    rlwe::SampleExtract(lwe_N_ct.data(), copy, extract_indices_,
                        runtime_->SEALRunTime());
  }
  lwe_ct.resize(parms_.nslots);

  ThreadPool pool(num_threads_);
  const size_t work_load = (parms_.nslots + num_threads_ - 1) / num_threads_;
  for (int w = 0; w < num_threads_; ++w) {
    size_t start = w * work_load;
    size_t end = std::min<size_t>(start + work_load, parms_.nslots);

    if (end > start) {
      pool.enqueue(
          [&](size_t s, size_t e) {
            for (size_t i = s; i < e; ++i) {
              rlwe::LWEKeySwitch(&lwe_ct[i], &lwe_N_ct[i], lvl2Tolvl0_,
                                 lvl0_runtime_);
              lwe_ct[i].scale = in.scale();
            }
          },
          start, end);
    }
  }
  return Status::Ok();
}

void PegasusRunTime::MulScalarLWECt(lwe::Ctx_st &out, const lwe::Ctx_st &a,
                                    uint64_t scalar) const {
  const auto &q0 =
      lvl0_runtime_->first_context_data()->parms().coeff_modulus()[0];
  using namespace seal::util;
  multiply_poly_scalar_coeffmod(CtData(&a), lwe::params::n() + 1, scalar, q0,
                                CtData(&out));
  out.scale = a.scale;
}

void PegasusRunTime::AddLWECt(lwe::Ctx_st &out, lwe::Ctx_st const &a,
                              lwe::Ctx_st const &b) const {
  using namespace lwe;
  using namespace seal::util;
  const auto &q0 =
      lvl0_runtime_->first_context_data()->parms().coeff_modulus()[0];
  if (!seal::util::are_close(a.scale, b.scale)) {
    throw std::invalid_argument("AddLWECt scale mismatch");
  }
  add_poly_coeffmod(CtData(&a), CtData(&b), lwe::params::n() + 1, q0,
                    CtData(&out));
  out.scale = a.scale;
}

void PegasusRunTime::SubLWECt(lwe::Ctx_st &out, lwe::Ctx_st const &a,
                              lwe::Ctx_st const &b) const {
  using namespace lwe;
  using namespace seal::util;
  if (&a == &b) {
    throw std::runtime_error("SubLWECt self-substraction");
  }
  if (!seal::util::are_close(a.scale, b.scale)) {
    throw std::invalid_argument("SubLWECt scale mismatch");
  }
  const auto &q0 =
      lvl0_runtime_->first_context_data()->parms().coeff_modulus()[0];
  sub_poly_coeffmod(CtData(&a), CtData(&b), lwe::params::n() + 1, q0,
                    CtData(&out));
  out.scale = a.scale;
}

void PegasusRunTime::MulConstant(lwe::Ctx_t lvl0_ct, double v) const {
  F64 out_scale = lvl0_ct->scale * lutFunctor_->GetPostMultiplier();
  rlwe::RLWE2LWECt_t lvl1_ct;
  lutFunctor_->MulConstant(lvl1_ct, lvl0_ct, v);
  rlwe::LWEKeySwitch(lvl0_ct, lvl1_ct, lvl1Tolvl0_, lvl0_runtime_);
  lvl0_ct->scale = out_scale;
}

double PegasusRunTime::DecryptLWE(lwe::Ctx_st const &lwe_ct) const {
  return lwe::SymDec(&lwe_ct, lvl0_sk_non_ntt_, lvl0_runtime_);
}

std::string PegasusRunTime::level2Params(int level_left) const {
  std::stringstream ss;
  ss << "{\"log2PolyDegree\":" << (int)std::log2(parms_.lvl2_lattice_dim)
     << ",";
  ss << "\"nSpecialPrimes\":1,\"seed\":0,";
  ss << "\"moduliArray\":[" << std::to_string(numBitsP0()) << ",";

  level_left = std::max(1, level_left);
  std::string norm_sze = std::to_string(numBitsP0());
  for (int i = 1; i < level_left; ++i) {
    ss << norm_sze << ",";
  }

  // Final is the special modulus
  ss << 59 << "]";
  ss << "}";
  return ss.str();
}

void PegasusRunTime::setUpRuntime() {
  using namespace seal;
  auto lvl2_runtime = runtime_->SEALRunTime();
  const auto &modulus =
      lvl2_runtime->key_context_data()->parms().coeff_modulus();
  EncryptionParameters parms(seal::scheme_type::CKKS);

  // Level-1 RGSW works with 2 moduli
  std::vector<Modulus> lvl1_modulus{modulus.front(), modulus.back()};
  parms.set_coeff_modulus(lvl1_modulus);
  parms.set_poly_modulus_degree(parms_.lvl1_lattice_dim);
  parms.set_galois_generator(5);
  lvl1_runtime_ = SEALContext::Create(parms, true, sec_level_type::none);

  // Level-0 LWE works with 1 modulus
  std::vector<Modulus> lvl0_modulus{modulus.front()};
  parms.set_poly_modulus_degree(parms_.lvl0_lattice_dim);
  parms.set_coeff_modulus(lvl0_modulus);
  lvl0_runtime_ = SEALContext::Create(parms, true, sec_level_type::none);
}

void PegasusRunTime::Binary(lwe::Ctx_t lvl0_ct, double threshold) const {
    F64 out_scale = lvl0_ct->scale * lutFunctor_->GetPostMultiplier();
    rlwe::RLWE2LWECt_t lvl1_ct;
    lutFunctor_->Binary(lvl1_ct, lvl0_ct, threshold);
    rlwe::LWEKeySwitch(lvl0_ct, lvl1_ct, lvl1Tolvl0_, lvl0_runtime_);
    lvl0_ct->scale = out_scale;
}

void PegasusRunTime::Binary(lwe::Ctx_st *lvl0_ct, size_t num_wires, double threshold) const {
    ThreadPool pool(num_threads_);
    const size_t work_load = (num_wires + num_threads_ - 1) / num_threads_;
    for (int w = 0; w < num_threads_; ++w) {
        size_t start = w * work_load;
        size_t end = std::min<size_t>(start + work_load, num_wires);
        if (end > start) {
            pool.enqueue(
                [&](size_t s, size_t e) {
                    for (size_t i = s; i < e; ++i) {
                        Binary(lvl0_ct + i, threshold);
                    }
                },
                start, end);
        }
    }
}

void PegasusRunTime::setUpSecretKeys() {
  // non-ntt form of the CKKS secret key
  auto lvl2_runtime = runtime_->SEALRunTime();
  auto const &lvl2_sk = runtime_->SEALSecretKey();
  lvl2_sk_non_ntt_.data().resize(parms_.lvl2_lattice_dim);
  std::copy_n((const uint64_t *)lvl2_sk.data().data(), parms_.lvl2_lattice_dim,
              lvl2_sk_non_ntt_.data().data());
  rlwe::SwitchNTTForm(lvl2_sk_non_ntt_.data().data(), NTTDir::FromNTT, 1,
                      lvl2_runtime);
  // Generate sk_{lut}
  lwe::GenerateHammingSecretKey(lvl1_sk_ntt_, SECRET_KEY_HW, /*is_ntt*/ true,
                                lvl1_runtime_);
  lvl1_sk_non_ntt_.data().resize(parms_.lvl1_lattice_dim);
  std::copy_n((const uint64_t *)lvl1_sk_ntt_.data().data(),
              parms_.lvl1_lattice_dim, lvl1_sk_non_ntt_.data().data());
  rlwe::SwitchNTTForm(lvl1_sk_non_ntt_.data().data(), NTTDir::FromNTT, 1,
                      lvl1_runtime_);
  // Generate sk_{lwe}
  lwe::SKInit(lvl0_sk_ntt_, lvl0_sk_non_ntt_, SECRET_KEY_HW, lvl0_runtime_);
}

void PegasusRunTime::setUpPublicKeys() {
  if (parms_.enable_repacking) {
    gemini::RpKeyInit(repackKey_, std::pow(2., numBitsP0()),
                      runtime_->SEALSecretKey(), lvl0_sk_non_ntt_, runtime_);
  }

  rlwe::BKInit(lutEvalKey_, lvl0_sk_non_ntt_, lvl1_sk_ntt_, lvl1_runtime_);
  rlwe::LWEKSKeyInit(lvl2Tolvl0_, KS_DC_BASE, lvl2_sk_non_ntt_.data(),
                     lvl0_sk_ntt_, lvl0_runtime_, runtime_->SEALRunTime());
  rlwe::LWEKSKeyInit(lvl1Tolvl0_, KS_DC_BASE, lvl1_sk_non_ntt_.data(),
                     lvl0_sk_ntt_, lvl0_runtime_, lvl1_runtime_);
}

void PegasusRunTime::setUpFunctors() {
  using namespace gemini;
  using namespace seal;
  LinearTransformer::Parms ltParams;
  ltParams.nslots = parms_.nslots;
  ltParams.s2c_lvl_start = 2;
  ltParams.c2s_lvl_start = 0;  // no C2S is needed
  ltParams.s2cMultiplier = parms_.s2c_multiplier;
  ltParams.c2sMultiplier = 1.;

  if (parms_.enable_repacking) {
    sinFunctor_.reset(new BetterSine(runtime_));
  }
  linearTransformer_.reset(new LinearTransformer(ltParams, runtime_));

  lutFunctor_.reset(
      new rlwe::LWEGateBooter(parms_.scale, lutEvalKey_, lvl1_runtime_, 1.));

  output_bounded_lutFunctor_.reset(new rlwe::LWEGateBooter(
      parms_.scale, lutEvalKey_, lvl1_runtime_, ExtraScaling()));

  extract_indices_.resize(parms_.nslots);
  const size_t log2N = (size_t)std::log2(parms_.lvl2_lattice_dim);
  for (size_t i = 0; i < parms_.nslots; ++i) {
    extract_indices_[i] = seal::util::reverse_bits(i, log2N - 1);
  }
}

Status PegasusRunTime::Repack(Ctx &out,
                              std::vector<lwe::Ctx_st> const &wires) const {
  if (!parms_.enable_repacking) {
    return Status::NotReady("Repacking is not ready");
  }

  size_t n = wires.size();
  if (n > parms_.lvl2_lattice_dim / 2) {
    return Status::ArgumentError("Repacking too many lwe ciphers to repack");
  }
  if (!IsTwoPower(n)) {
    return Status::ArgumentError("Repacking invalid number of lwe ciphers");
  }

  CHECK_STATUS(::gemini::Repack(out, sinFunctor_->interval_bound(), wires,
                                repackKey_, runtime_));
  out.scale() = parms_.scale;
  CHECK_STATUS(sinFunctor_->Apply(out, parms_.scale));
  return Status::Ok();
}

bool PegasusRunTime::SaveSealContexts(std::ostream &os) const {

    if (runtime_) {
        bool has_runtime = true;
        os.write(reinterpret_cast<const char*>(&has_runtime), sizeof(has_runtime));
        if (!runtime_->Save(os)) {
            return false;
        }
    } else {
        bool has_runtime = false;
        os.write(reinterpret_cast<const char*>(&has_runtime), sizeof(has_runtime));
    }

    return true;
}

bool PegasusRunTime::LoadSealContexts(std::istream &is, std::string const& json) {
    bool has_runtime;
    
    is.read(reinterpret_cast<char*>(&has_runtime), sizeof(has_runtime));
    if (!is.good()) return false;

    if (has_runtime) {
        if (!runtime_->Load(is, json)) {
            return false;
        }
      std::cout << "Load runtime success" << std::endl;
    } else {
        runtime_ = nullptr;
    }

    return true;
}

Status PegasusRunTime::Save(std::ostream &os) const {
    int version = 1;
    os.write(reinterpret_cast<const char*>(&version), sizeof(version));

    os.write(reinterpret_cast<const char*>(&parms_), sizeof(parms_));
    if (!os.good()) return Status::FileIOError("Error writing Parms");

    os.write(reinterpret_cast<const char*>(&num_threads_), sizeof(num_threads_));
    if (!os.good()) return Status::FileIOError("Error writing num_threads_");

    if (!SaveSealContexts(os)) {
        return Status::FileIOError("Error saving SEAL contexts and RunTime");
    }

    if (!os.good()) {
        return Status::FileIOError("Error occurred during saving PegasusRunTime");
    }

    return Status::Ok();
}

Status PegasusRunTime::Load(std::istream &is, std::string const& json) {
    if (!is.good()) {
        return Status::FileIOError("Input stream not good for reading");
    }
    int version;
    is.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (!is.good()) return Status::FileIOError("Error reading version");
    if (version != 1) {
        return Status::ArgumentError("Unsupported version");
    }
    is.read(reinterpret_cast<char*>(&parms_), sizeof(parms_));
    if (!is.good()) return Status::FileIOError("Error reading Parms");

    size_t &non_const_num_threads = const_cast<size_t&>(num_threads_);
    is.read(reinterpret_cast<char*>(&non_const_num_threads), sizeof(non_const_num_threads));
    if (!is.good()) return Status::FileIOError("Error reading num_threads_");
    if (!LoadSealContexts(is, json)) {
        return Status::FileIOError("Error loading SEAL contexts and RunTime");
    }

    if (!is.good()) {
        return Status::FileIOError("Error occurred during loading PegasusRunTime");
    }

    return Status::Ok();
}

Status PegasusRunTime::SaveToFile(const std::string &filename) const {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs.is_open()) {
        return Status::FileIOError("Failed to open file for writing");
    }

    return Save(ofs);
}

Status PegasusRunTime::LoadFromFile(const std::string &filename, std::string const& json) {
    std::ifstream ifs(filename, std::ios::binary);
    if (!ifs.is_open()) {
        return Status::FileIOError("Failed to open file for reading");
    }
    Status load_status = Load(ifs, json);
    std::cout << "LoadFromFile status: " << load_status << std::endl;
    return Status::Ok();
}


}  // namespace gemini

