#include "pegasus/pegasus_runtime.h"
#include "pegasus/timer.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <type_traits>
using namespace gemini;
using namespace std;


Ctx rlwe_multiply(Ctx &a, Ctx &b, PegasusRunTime &pg_rt)
{
    Ctx result = a;
    pg_rt.Mul(result, b);
    pg_rt.RelinThenRescale(result);
    return result;
}

Ctx rlwe_addition(Ctx &a, Ctx &b, PegasusRunTime &pg_rt)
{
    Ctx result = a;
    pg_rt.Add(result, b);
    return result;
}

Ctx rlwe_substraction(Ctx &a, Ctx &b, PegasusRunTime &pg_rt)
{
    Ctx result = a;
    pg_rt.Sub(result, b);
    return result;
}

Ctx RotateLeft(Ctx &a, int step, PegasusRunTime &pg_rt)
{
    Ctx result = a;
    pg_rt.RotateLeft(result, (size_t)abs(step));
    return result;
}
        
std::vector<lwe::Ctx_st> binary64x64(Ctx &v1, std::vector<lwe::Ctx_st> &v2, PegasusRunTime &pg_rt)
{
    std::vector<lwe::Ctx_st> v3;
    pg_rt.ExtraAllCoefficients(v1, v3);
    pg_rt.Binary(v3.data(), v3.size(), 0.49019607843137253);
    v2 = v3;
    return v2;
}

int main() {
    PegasusRunTime::Parms pp;
    pp.lvl0_lattice_dim = lwe::params::n();
    pp.lvl1_lattice_dim = 1 << 12;
    pp.lvl2_lattice_dim = 1 << 16;
    pp.nlevels = 4;
    pp.scale = std::pow(2., 40);
    pp.nslots = 64;
    pp.s2c_multiplier = 1.;
    pp.enable_repacking = false;

    // PegasusRunTime pg_rt(pp, 4);
    std::string filename = "/home/shenruiyu/fhetran/fhe-transpiler-demo/thirdparty/OpenPEGASUS/examples/runtime.bin";
    std::string filename2 = "/home/shenruiyu/fhetran/fhe-transpiler-demo/thirdparty/OpenPEGASUS/examples/ciphertext.bin";
    // gemini::Status save_status = pg_rt.SaveToFile(filename);
    // std::cout << "PegasusRunTime saved successfully to " << filename << std::endl;
    // F64Vec slots(64, 0.8);
    // slots[0] = 0.0;
    // Ctx v1;
    // pg_rt.EncodeThenEncrypt(slots, v1);
    // std::ofstream ofs(filename2, std::ios::binary);
    // v1.save(ofs);
    // ofs.close();


    PegasusRunTime pgrt_loaded(pp, 4, filename);
    std::cout << "PegasusRunTime loaded successfully from " << filename << std::endl;

    Ctx v10;
    std::ifstream ifs(filename2, std::ios::binary);
    pgrt_loaded.load_cipher(v10, ifs);
    ifs.close();

    F64Vec slots3;
    




    auto v4 = rlwe_addition(v10, v10, pgrt_loaded);
    auto v5 = RotateLeft(v4, 1, pgrt_loaded);
    pgrt_loaded.DecryptThenDecode(v5, slots3);
    for (int i = 0; i < 64; i++) {
        std::cout << slots3[i] << " ";
    }


    // std::vector<lwe::Ctx_st> v2;
    // pgrt_loaded.SlotsToCoeffs(v10);
    // std::vector<lwe::Ctx_st> result = binary64x64(v10, v2, pgrt_loaded);

    // F64Vec output;
    // for (auto &lwe_ct : result) {
    //     double value = pgrt_loaded.DecryptLWE(lwe_ct);
    //     output.push_back(value);
    // }

    // std::ofstream outfile("/home/shenruiyu/fhetran/fhe-transpiler-demo/test/binary64x64/output_image.txt");
    // if (outfile.is_open()) {
    //     size_t index = 0;
    //     for (size_t i = 0; i < 8; ++i) {
    //         for (size_t j = 0; j < 8; ++j) {
    //             outfile << output[index++]*255 << " ";
    //         }
    //         outfile << std::endl;
    //     }
    //     outfile.close();
    //     std::cout << "Data saved in output_image.txt" << std::endl;
    // } else {
    //     std::cerr << "Can't write" << std::endl;
    // }





    // F64Vec slots2(64, 0.0);
    // pgrt_loaded.DecryptThenDecode(v5, slots2);
    // for (int i = 0; i < 64; i++) {
    //     std::cout << slots2[i] << " ";
    // }


    return 0;
}