#include "chem/ligand_factory.hpp"

#include <iterator>

namespace scheme {
namespace chemical {

struct AtomIncludeCheck {
  bool hydrogen;
  AtomIncludeCheck(bool h) : hydrogen(h) {}

  bool operator()(int at, int nc, int bbc) const {
    if (at == 21) return false;
    if (!hydrogen && at > 23) return false;
    if (bbc < nc) return false;
    if (bbc > nc && at > 23) return false;
    return true;
  }
};

std::vector<std::string> get_pdb_lines(std::string resn, bool hydrogen) {
  using std::vector;
  using std::string;

  vector<std::string> lines;
  std::back_insert_iterator<vector<string>> outiter(lines);

  AtomIncludeCheck check(hydrogen);

  if (resn == "ALA") {
    vector<int> NC(10);
    vector<int> BBC(10);
    vector<int> AT(10); /*ALA_ N  */
    AT[0] = 18;
    NC[0] = 0;
    BBC[0] = 0;
    /*ALA_ CA */ AT[1] = 19;
    NC[1] = 0;
    BBC[1] = 0;
    /*ALA_ C  */ AT[2] = 20;
    NC[2] = 0;
    BBC[2] = 0;
    /*ALA_ O  */ AT[3] = 21;
    NC[3] = 0;
    BBC[3] = 0;
    /*ALA_ CB */ AT[4] = 5;
    NC[4] = 0;
    BBC[4] = 0;
    /*ALA_ H  */ AT[5] = 27;
    NC[5] = 0;
    BBC[5] = 0;
    /*ALA_ HA */ AT[6] = 25;
    NC[6] = 0;
    BBC[6] = 0;
    /*ALA_1HB */ AT[7] = 25;
    NC[7] = 0;
    BBC[7] = 0;
    /*ALA_2HB */ AT[8] = 25;
    NC[8] = 0;
    BBC[8] = 0;
    /*ALA_3HB */ AT[9] = 25;
    NC[9] = 0;
    BBC[9] = 0;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   ALA     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  ALA     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   ALA     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   ALA     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  ALA     0       1.988  -0.773  -1.199");
    ++j;
    // if( check(AT[j],NC[j],BBC[j]) ) *outiter++ = std::string("ATOM      6  H
    // ALA     0      -0.492   0.764  -0.441"); ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  HA  ALA     0       1.797  -0.490   0.913");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8 1HB  ALA     0       3.078  -0.764  -1.185");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9 2HB  ALA     0       1.633  -1.802  -1.154");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10 3HB  ALA     0       1.633  -0.307  -2.117");
    ++j;
  }
  if (resn == "CYS") {
    vector<int> NC(11);
    vector<int> BBC(11);
    vector<int> AT(11); /*CYS_ N  */
    AT[0] = 18;
    NC[0] = 1;
    BBC[0] = 0;
    /*CYS_ CA */ AT[1] = 19;
    NC[1] = 1;
    BBC[1] = 2;
    /*CYS_ C  */ AT[2] = 20;
    NC[2] = 1;
    BBC[2] = 0;
    /*CYS_ O  */ AT[3] = 21;
    NC[3] = 1;
    BBC[3] = 0;
    /*CYS_ CB */ AT[4] = 4;
    NC[4] = 1;
    BBC[4] = 2;
    /*CYS_ SG */ AT[5] = 17;
    NC[5] = 1;
    BBC[5] = 2;
    /*CYS_ H  */ AT[6] = 27;
    NC[6] = 1;
    BBC[6] = 0;
    /*CYS_ HA */ AT[7] = 25;
    NC[7] = 1;
    BBC[7] = 0;
    /*CYS_1HB */ AT[8] = 25;
    NC[8] = 1;
    BBC[8] = 1;
    /*CYS_2HB */ AT[9] = 25;
    NC[9] = 1;
    BBC[9] = 1;
    /*CYS_ HG */ AT[10] = 25;
    NC[10] = 1;
    BBC[10] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   CYS     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  CYS     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   CYS     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   CYS     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  CYS     0       1.996  -0.750  -1.219");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  SG  CYS     0       0.710  -1.417  -2.303");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  H   CYS     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  HA  CYS     0       1.804  -0.507   0.901");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9 1HB  CYS     0       2.622  -0.082  -1.810");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10 2HB  CYS     0       2.622  -1.579  -0.889");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  HG  CYS     0       1.542  -1.960  -3.186");
    ++j;
  }
  if (resn == "ASP") {
    vector<int> NC(12);
    vector<int> BBC(12);
    vector<int> AT(12); /*ASP_ N  */
    AT[0] = 18;
    NC[0] = 2;
    BBC[0] = 0;
    /*ASP_ CA */ AT[1] = 19;
    NC[1] = 2;
    BBC[1] = 0;
    /*ASP_ C  */ AT[2] = 20;
    NC[2] = 2;
    BBC[2] = 0;
    /*ASP_ O  */ AT[3] = 21;
    NC[3] = 2;
    BBC[3] = 0;
    /*ASP_ CB */ AT[4] = 4;
    NC[4] = 2;
    BBC[4] = 2;
    /*ASP_ CG */ AT[5] = 2;
    NC[5] = 2;
    BBC[5] = 2;
    /*ASP_ OD1*/ AT[6] = 15;
    NC[6] = 2;
    BBC[6] = 2;
    /*ASP_ OD2*/ AT[7] = 15;
    NC[7] = 2;
    BBC[7] = 2;
    /*ASP_ H  */ AT[8] = 27;
    NC[8] = 2;
    BBC[8] = 0;
    /*ASP_ HA */ AT[9] = 25;
    NC[9] = 2;
    BBC[9] = 0;
    /*ASP_1HB */ AT[10] = 25;
    NC[10] = 2;
    BBC[10] = 1;
    /*ASP_2HB */ AT[11] = 25;
    NC[11] = 2;
    BBC[11] = 1;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   ASP     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  ASP     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   ASP     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   ASP     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  ASP     0       1.995  -0.762  -1.214");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  ASP     0       0.888  -1.318  -2.101");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  OD1 ASP     0      -0.259  -1.118  -1.782");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  OD2 ASP     0       1.202  -1.938  -3.089");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  H   ASP     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  HA  ASP     0       1.804  -0.500   0.905");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11 1HB  ASP     0       2.621  -0.099  -1.812");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 2HB  ASP     0       2.621  -1.589  -0.878");
    ++j;
  }
  if (resn == "GLU") {
    vector<int> NC(15);
    vector<int> BBC(15);
    vector<int> AT(15); /*GLU_ N  */
    AT[0] = 18;
    NC[0] = 3;
    BBC[0] = 0;
    /*GLU_ CA */ AT[1] = 19;
    NC[1] = 3;
    BBC[1] = 0;
    /*GLU_ C  */ AT[2] = 20;
    NC[2] = 3;
    BBC[2] = 0;
    /*GLU_ O  */ AT[3] = 21;
    NC[3] = 3;
    BBC[3] = 0;
    /*GLU_ CB */ AT[4] = 4;
    NC[4] = 3;
    BBC[4] = 2;
    /*GLU_ CG */ AT[5] = 4;
    NC[5] = 3;
    BBC[5] = 3;
    /*GLU_ CD */ AT[6] = 2;
    NC[6] = 3;
    BBC[6] = 3;
    /*GLU_ OE1*/ AT[7] = 15;
    NC[7] = 3;
    BBC[7] = 3;
    /*GLU_ OE2*/ AT[8] = 15;
    NC[8] = 3;
    BBC[8] = 3;
    /*GLU_ H  */ AT[9] = 27;
    NC[9] = 3;
    BBC[9] = 0;
    /*GLU_ HA */ AT[10] = 25;
    NC[10] = 3;
    BBC[10] = 0;
    /*GLU_1HB */ AT[11] = 25;
    NC[11] = 3;
    BBC[11] = 1;
    /*GLU_2HB */ AT[12] = 25;
    NC[12] = 3;
    BBC[12] = 1;
    /*GLU_1HG */ AT[13] = 25;
    NC[13] = 3;
    BBC[13] = 2;
    /*GLU_2HG */ AT[14] = 25;
    NC[14] = 3;
    BBC[14] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   GLU     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  GLU     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   GLU     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   GLU     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  GLU     0       1.991  -0.764  -1.214");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  GLU     0       0.911  -1.336  -2.121");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CD  GLU     0      -0.480  -1.032  -1.638");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  OE1 GLU     0      -0.609  -0.392  -0.622");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  OE2 GLU     0      -1.414  -1.440  -2.287");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  H   GLU     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  HA  GLU     0       1.804  -0.498   0.906");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 1HB  GLU     0       2.617  -0.103  -1.814");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 2HB  GLU     0       2.617  -1.591  -0.877");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 1HG  GLU     0       1.035  -0.921  -3.121");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 2HG  GLU     0       1.040  -2.415  -2.187");
    ++j;
  }
  if (resn == "PHE") {
    vector<int> NC(20);
    vector<int> BBC(20);
    vector<int> AT(20); /*PHE_ N  */
    AT[0] = 18;
    NC[0] = 2;
    BBC[0] = 0;
    /*PHE_ CA */ AT[1] = 19;
    NC[1] = 2;
    BBC[1] = 0;
    /*PHE_ C  */ AT[2] = 20;
    NC[2] = 2;
    BBC[2] = 0;
    /*PHE_ O  */ AT[3] = 21;
    NC[3] = 2;
    BBC[3] = 0;
    /*PHE_ CB */ AT[4] = 4;
    NC[4] = 2;
    BBC[4] = 2;
    /*PHE_ CG */ AT[5] = 6;
    NC[5] = 2;
    BBC[5] = 2;
    /*PHE_ CD1*/ AT[6] = 6;
    NC[6] = 2;
    BBC[6] = 2;
    /*PHE_ CD2*/ AT[7] = 6;
    NC[7] = 2;
    BBC[7] = 2;
    /*PHE_ CE1*/ AT[8] = 6;
    NC[8] = 2;
    BBC[8] = 2;
    /*PHE_ CE2*/ AT[9] = 6;
    NC[9] = 2;
    BBC[9] = 2;
    /*PHE_ CZ */ AT[10] = 6;
    NC[10] = 2;
    BBC[10] = 2;
    /*PHE_ H  */ AT[11] = 27;
    NC[11] = 2;
    BBC[11] = 0;
    /*PHE_ HA */ AT[12] = 25;
    NC[12] = 2;
    BBC[12] = 0;
    /*PHE_1HB */ AT[13] = 25;
    NC[13] = 2;
    BBC[13] = 1;
    /*PHE_2HB */ AT[14] = 25;
    NC[14] = 2;
    BBC[14] = 1;
    /*PHE_ HD1*/ AT[15] = 26;
    NC[15] = 2;
    BBC[15] = 2;
    /*PHE_ HD2*/ AT[16] = 26;
    NC[16] = 2;
    BBC[16] = 2;
    /*PHE_ HE1*/ AT[17] = 26;
    NC[17] = 2;
    BBC[17] = 2;
    /*PHE_ HE2*/ AT[18] = 26;
    NC[18] = 2;
    BBC[18] = 2;
    /*PHE_ HZ */ AT[19] = 26;
    NC[19] = 2;
    BBC[19] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   PHE     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  PHE     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   PHE     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   PHE     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  PHE     0       1.992  -0.760  -1.216");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  PHE     0       0.916  -1.315  -2.104");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CD1 PHE     0      -0.423  -1.124  -1.798");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  CD2 PHE     0       1.240  -2.030  -3.247");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  CE1 PHE     0      -1.414  -1.634  -2.615");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  CE2 PHE     0       0.252  -2.540  -4.066");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  CZ  PHE     0      -1.077  -2.342  -3.749");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12  H   PHE     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13  HA  PHE     0       1.804  -0.502   0.905");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 1HB  PHE     0       2.618  -0.097  -1.812");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 2HB  PHE     0       2.618  -1.586  -0.882");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16  HD1 PHE     0      -0.690  -0.563  -0.902");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17  HD2 PHE     0       2.290  -2.187  -3.497");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     18  HE1 PHE     0      -2.463  -1.477  -2.363");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     19  HE2 PHE     0       0.520  -3.100  -4.962");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     20  HZ  PHE     0      -1.857  -2.746  -4.393");
    ++j;
  }
  if (resn == "GLY") {
    vector<int> NC(7);
    vector<int> BBC(7);
    vector<int> AT(7); /*GLY_ N  */
    AT[0] = 18;
    NC[0] = 0;
    BBC[0] = 0;
    /*GLY_ CA */ AT[1] = 19;
    NC[1] = 0;
    BBC[1] = 0;
    /*GLY_ C  */ AT[2] = 20;
    NC[2] = 0;
    BBC[2] = 0;
    /*GLY_ O  */ AT[3] = 21;
    NC[3] = 0;
    BBC[3] = 0;
    /*GLY_ H  */ AT[4] = 27;
    NC[4] = 0;
    BBC[4] = 0;
    /*GLY_1HA */ AT[5] = 25;
    NC[5] = 0;
    BBC[5] = 0;
    /*GLY_2HA */ AT[6] = 25;
    NC[6] = 0;
    BBC[6] = 0;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   GLY     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  GLY     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   GLY     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   GLY     0       1.383   2.339  -0.529");
    ++j;
    // if( check(AT[j],NC[j],BBC[j]) ) *outiter++ = std::string("ATOM      5  H
    // GLY     0      -0.492   0.764  -0.441"); ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6 1HA  GLY     0       1.822  -0.535   0.877");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7 2HA  GLY     0       1.822  -0.535  -0.876");
    ++j;
  }
  if (resn == "HIS") {
    vector<int> NC(17);
    vector<int> BBC(17);
    vector<int> AT(17); /*HIS_ N  */
    AT[0] = 18;
    NC[0] = 2;
    BBC[0] = 0;
    /*HIS_ CA */ AT[1] = 19;
    NC[1] = 2;
    BBC[1] = 0;
    /*HIS_ C  */ AT[2] = 20;
    NC[2] = 2;
    BBC[2] = 0;
    /*HIS_ O  */ AT[3] = 21;
    NC[3] = 2;
    BBC[3] = 0;
    /*HIS_ CB */ AT[4] = 4;
    NC[4] = 2;
    BBC[4] = 2;
    /*HIS_ CG */ AT[5] = 6;
    NC[5] = 2;
    BBC[5] = 2;
    /*HIS_ ND1*/ AT[6] = 8;
    NC[6] = 2;
    BBC[6] = 2;
    /*HIS_ CD2*/ AT[7] = 6;
    NC[7] = 2;
    BBC[7] = 2;
    /*HIS_ CE1*/ AT[8] = 6;
    NC[8] = 2;
    BBC[8] = 2;
    /*HIS_ NE2*/ AT[9] = 7;
    NC[9] = 2;
    BBC[9] = 2;
    /*HIS_ H  */ AT[10] = 27;
    NC[10] = 2;
    BBC[10] = 0;
    /*HIS_ HA */ AT[11] = 25;
    NC[11] = 2;
    BBC[11] = 0;
    /*HIS_1HB */ AT[12] = 25;
    NC[12] = 2;
    BBC[12] = 1;
    /*HIS_2HB */ AT[13] = 25;
    NC[13] = 2;
    BBC[13] = 1;
    /*HIS_ HD2*/ AT[14] = 25;
    NC[14] = 2;
    BBC[14] = 2;
    /*HIS_ HE1*/ AT[15] = 25;
    NC[15] = 2;
    BBC[15] = 2;
    /*HIS_ HE2*/ AT[16] = 24;
    NC[16] = 2;
    BBC[16] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   HIS     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  HIS     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   HIS     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   HIS     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  HIS     0       1.999  -0.764  -1.213");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  HIS     0       0.929  -1.322  -2.099");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  ND1 HIS     0      -0.415  -1.157  -1.837");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  CD2 HIS     0       1.004  -2.043  -3.242");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  CE1 HIS     0      -1.121  -1.752  -2.782");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  NE2 HIS     0      -0.283  -2.297  -3.646");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  H   HIS     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12  HA  HIS     0       1.822  -0.495   0.900");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 1HB  HIS     0       2.626  -0.101  -1.810");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 2HB  HIS     0       2.626  -1.588  -0.874");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15  HD2 HIS     0       1.916  -2.362  -3.748");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16  HE1 HIS     0      -2.208  -1.787  -2.840");
    ++j;
    // if( check(AT[j],NC[j],BBC[j]) ) *outiter++ = std::string("ATOM     17
    // HE2 HIS     0      -0.542  -2.818  -4.472"); ++j;
  }
  if (resn == "ILE") {
    vector<int> NC(19);
    vector<int> BBC(19);
    vector<int> AT(19); /*ILE_ N  */
    AT[0] = 18;
    NC[0] = 2;
    BBC[0] = 0;
    /*ILE_ CA */ AT[1] = 19;
    NC[1] = 2;
    BBC[1] = 0;
    /*ILE_ C  */ AT[2] = 20;
    NC[2] = 2;
    BBC[2] = 0;
    /*ILE_ O  */ AT[3] = 21;
    NC[3] = 2;
    BBC[3] = 0;
    /*ILE_ CB */ AT[4] = 3;
    NC[4] = 2;
    BBC[4] = 2;
    /*ILE_ CG1*/ AT[5] = 4;
    NC[5] = 2;
    BBC[5] = 2;
    /*ILE_ CG2*/ AT[6] = 5;
    NC[6] = 2;
    BBC[6] = 1;
    /*ILE_ CD1*/ AT[7] = 5;
    NC[7] = 2;
    BBC[7] = 2;
    /*ILE_ H  */ AT[8] = 27;
    NC[8] = 2;
    BBC[8] = 0;
    /*ILE_ HA */ AT[9] = 25;
    NC[9] = 2;
    BBC[9] = 0;
    /*ILE_ HB */ AT[10] = 25;
    NC[10] = 2;
    BBC[10] = 1;
    /*ILE_1HG1*/ AT[11] = 25;
    NC[11] = 2;
    BBC[11] = 2;
    /*ILE_2HG1*/ AT[12] = 25;
    NC[12] = 2;
    BBC[12] = 2;
    /*ILE_1HG2*/ AT[13] = 25;
    NC[13] = 2;
    BBC[13] = 1;
    /*ILE_2HG2*/ AT[14] = 25;
    NC[14] = 2;
    BBC[14] = 1;
    /*ILE_3HG2*/ AT[15] = 25;
    NC[15] = 2;
    BBC[15] = 1;
    /*ILE_1HD1*/ AT[16] = 25;
    NC[16] = 2;
    BBC[16] = 2;
    /*ILE_2HD1*/ AT[17] = 25;
    NC[17] = 2;
    BBC[17] = 2;
    /*ILE_3HD1*/ AT[18] = 25;
    NC[18] = 2;
    BBC[18] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   ILE     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  ILE     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   ILE     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   ILE     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  ILE     0       2.007  -0.764  -1.218");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG1 ILE     0       0.857  -1.301  -2.075");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CG2 ILE     0       2.916  -1.899  -0.770");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  CD1 ILE     0      -0.514  -0.964  -1.536");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  H   ILE     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  HA  ILE     0       1.804  -0.499   0.904");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  HB  ILE     0       2.577  -0.083  -1.849");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 1HG1 ILE     0       0.937  -0.899  -3.084");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 2HG1 ILE     0       0.935  -2.386  -2.149");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 1HG2 ILE     0       3.295  -2.429  -1.643");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 2HG2 ILE     0       3.751  -1.492  -0.201");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16 3HG2 ILE     0       2.352  -2.589  -0.142");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17 1HD1 ILE     0      -1.277  -1.378  -2.196");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     18 2HD1 ILE     0      -0.628  -1.388  -0.538");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     19 3HD1 ILE     0      -0.629   0.118  -1.486");
    ++j;
  }
  if (resn == "LYS") {
    vector<int> NC(22);
    vector<int> BBC(22);
    vector<int> AT(22); /*LYS_ N  */
    AT[0] = 18;
    NC[0] = 4;
    BBC[0] = 0;
    /*LYS_ CA */ AT[1] = 19;
    NC[1] = 4;
    BBC[1] = 0;
    /*LYS_ C  */ AT[2] = 20;
    NC[2] = 4;
    BBC[2] = 2;
    /*LYS_ O  */ AT[3] = 21;
    NC[3] = 4;
    BBC[3] = 2;
    /*LYS_ CB */ AT[4] = 4;
    NC[4] = 4;
    BBC[4] = 2;
    /*LYS_ CG */ AT[5] = 4;
    NC[5] = 4;
    BBC[5] = 3;
    /*LYS_ CD */ AT[6] = 4;
    NC[6] = 4;
    BBC[6] = 4;
    /*LYS_ CE */ AT[7] = 4;
    NC[7] = 4;
    BBC[7] = 4;
    /*LYS_ NZ */ AT[8] = 10;
    NC[8] = 4;
    BBC[8] = 4;
    /*LYS_ H  */ AT[9] = 27;
    NC[9] = 4;
    BBC[9] = 2;
    /*LYS_ HA */ AT[10] = 25;
    NC[10] = 4;
    BBC[10] = 2;
    /*LYS_1HB */ AT[11] = 25;
    NC[11] = 4;
    BBC[11] = 2;
    /*LYS_2HB */ AT[12] = 25;
    NC[12] = 4;
    BBC[12] = 2;
    /*LYS_1HG */ AT[13] = 25;
    NC[13] = 4;
    BBC[13] = 2;
    /*LYS_2HG */ AT[14] = 25;
    NC[14] = 4;
    BBC[14] = 2;
    /*LYS_1HD */ AT[15] = 25;
    NC[15] = 4;
    BBC[15] = 3;
    /*LYS_2HD */ AT[16] = 25;
    NC[16] = 4;
    BBC[16] = 3;
    /*LYS_1HE */ AT[17] = 25;
    NC[17] = 4;
    BBC[17] = 4;
    /*LYS_2HE */ AT[18] = 25;
    NC[18] = 4;
    BBC[18] = 4;
    /*LYS_1HZ */ AT[19] = 24;
    NC[19] = 4;
    BBC[19] = 4;
    /*LYS_2HZ */ AT[20] = 24;
    NC[20] = 4;
    BBC[20] = 4;
    /*LYS_3HZ */ AT[21] = 24;
    NC[21] = 4;
    BBC[21] = 4;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   LYS     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  LYS     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   LYS     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   LYS     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  LYS     0       1.994  -0.772  -1.207");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  LYS     0       0.915  -1.351  -2.113");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CD  LYS     0      -0.477  -1.020  -1.595");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  CE  LYS     0      -0.412  -0.201  -0.315");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  NZ  LYS     0       0.990   0.067   0.105");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  H   LYS     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  HA  LYS     0       1.804  -0.492   0.910");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 1HB  LYS     0       2.620  -0.114  -1.811");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 2HB  LYS     0       2.620  -1.596  -0.863");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 1HG  LYS     0       1.026  -0.943  -3.118");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 2HG  LYS     0       1.026  -2.434  -2.164");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16 1HD  LYS     0      -1.022  -0.454  -2.352");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17 2HD  LYS     0      -1.021  -1.944  -1.397");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     18 1HE  LYS     0      -0.922   0.750  -0.465");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     19 2HE  LYS     0      -0.921  -0.738   0.486");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     20 1HZ  LYS     0       0.990   0.611   0.956");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     21 2HZ  LYS     0       1.467  -0.809   0.265");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     22 3HZ  LYS     0       1.467   0.581  -0.622");
    ++j;
  }
  if (resn == "LEU") {
    vector<int> NC(19);
    vector<int> BBC(19);
    vector<int> AT(19); /*LEU_ N  */
    AT[0] = 18;
    NC[0] = 2;
    BBC[0] = 0;
    /*LEU_ CA */ AT[1] = 19;
    NC[1] = 2;
    BBC[1] = 0;
    /*LEU_ C  */ AT[2] = 20;
    NC[2] = 2;
    BBC[2] = 0;
    /*LEU_ O  */ AT[3] = 21;
    NC[3] = 2;
    BBC[3] = 0;
    /*LEU_ CB */ AT[4] = 4;
    NC[4] = 2;
    BBC[4] = 2;
    /*LEU_ CG */ AT[5] = 3;
    NC[5] = 2;
    BBC[5] = 2;
    /*LEU_ CD1*/ AT[6] = 5;
    NC[6] = 2;
    BBC[6] = 2;
    /*LEU_ CD2*/ AT[7] = 5;
    NC[7] = 2;
    BBC[7] = 2;
    /*LEU_ H  */ AT[8] = 27;
    NC[8] = 2;
    BBC[8] = 0;
    /*LEU_ HA */ AT[9] = 25;
    NC[9] = 2;
    BBC[9] = 0;
    /*LEU_1HB */ AT[10] = 25;
    NC[10] = 2;
    BBC[10] = 1;
    /*LEU_2HB */ AT[11] = 25;
    NC[11] = 2;
    BBC[11] = 1;
    /*LEU_ HG */ AT[12] = 25;
    NC[12] = 2;
    BBC[12] = 2;
    /*LEU_1HD1*/ AT[13] = 25;
    NC[13] = 2;
    BBC[13] = 2;
    /*LEU_2HD1*/ AT[14] = 25;
    NC[14] = 2;
    BBC[14] = 2;
    /*LEU_3HD1*/ AT[15] = 25;
    NC[15] = 2;
    BBC[15] = 2;
    /*LEU_1HD2*/ AT[16] = 25;
    NC[16] = 2;
    BBC[16] = 2;
    /*LEU_2HD2*/ AT[17] = 25;
    NC[17] = 2;
    BBC[17] = 2;
    /*LEU_3HD2*/ AT[18] = 25;
    NC[18] = 2;
    BBC[18] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   LEU     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  LEU     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   LEU     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   LEU     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  LEU     0       1.988  -0.761  -1.222");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  LEU     0       0.921  -1.343  -2.157");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CD1 LEU     0      -0.464  -1.008  -1.620");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  CD2 LEU     0       1.113  -0.783  -3.559");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  H   LEU     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  HA  LEU     0       1.804  -0.502   0.903");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11 1HB  LEU     0       2.610  -0.086  -1.808");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 2HB  LEU     0       2.610  -1.585  -0.875");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13  HG  LEU     0       1.014  -2.429  -2.185");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 1HD1 LEU     0      -1.222  -1.422  -2.284");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 2HD1 LEU     0      -0.581  -1.437  -0.625");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16 3HD1 LEU     0      -0.581   0.074  -1.565");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17 1HD2 LEU     0       0.355  -1.197  -4.224");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     18 2HD2 LEU     0       1.019   0.303  -3.533");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     19 3HD2 LEU     0       2.104  -1.052  -3.925");
    ++j;
  }
  if (resn == "MET") {
    vector<int> NC(17);
    vector<int> BBC(17);
    vector<int> AT(17); /*MET_ N  */
    AT[0] = 18;
    NC[0] = 3;
    BBC[0] = 0;
    /*MET_ CA */ AT[1] = 19;
    NC[1] = 3;
    BBC[1] = 0;
    /*MET_ C  */ AT[2] = 20;
    NC[2] = 3;
    BBC[2] = 0;
    /*MET_ O  */ AT[3] = 21;
    NC[3] = 3;
    BBC[3] = 0;
    /*MET_ CB */ AT[4] = 4;
    NC[4] = 3;
    BBC[4] = 2;
    /*MET_ CG */ AT[5] = 4;
    NC[5] = 3;
    BBC[5] = 3;
    /*MET_ SD */ AT[6] = 17;
    NC[6] = 3;
    BBC[6] = 3;
    /*MET_ CE */ AT[7] = 5;
    NC[7] = 3;
    BBC[7] = 3;
    /*MET_ H  */ AT[8] = 27;
    NC[8] = 3;
    BBC[8] = 0;
    /*MET_ HA */ AT[9] = 25;
    NC[9] = 3;
    BBC[9] = 0;
    /*MET_1HB */ AT[10] = 25;
    NC[10] = 3;
    BBC[10] = 1;
    /*MET_2HB */ AT[11] = 25;
    NC[11] = 3;
    BBC[11] = 1;
    /*MET_1HG */ AT[12] = 25;
    NC[12] = 3;
    BBC[12] = 2;
    /*MET_2HG */ AT[13] = 25;
    NC[13] = 3;
    BBC[13] = 2;
    /*MET_1HE */ AT[14] = 25;
    NC[14] = 3;
    BBC[14] = 3;
    /*MET_2HE */ AT[15] = 25;
    NC[15] = 3;
    BBC[15] = 3;
    /*MET_3HE */ AT[16] = 25;
    NC[16] = 3;
    BBC[16] = 3;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   MET     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  MET     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   MET     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   MET     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  MET     0       1.986  -0.776  -1.205");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  MET     0       0.903  -1.356  -2.104");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  SD  MET     0      -0.762  -0.979  -1.519");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  CE  MET     0      -0.396  -0.030  -0.046");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  H   MET     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  HA  MET     0       1.804  -0.488   0.911");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11 1HB  MET     0       2.610  -0.124  -1.813");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 2HB  MET     0       2.611  -1.601  -0.861");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 1HG  MET     0       1.013  -0.955  -3.111");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 2HG  MET     0       1.013  -2.439  -2.154");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 1HE  MET     0      -1.328   0.277   0.430");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16 2HE  MET     0       0.182  -0.642   0.647");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17 3HE  MET     0       0.181   0.855  -0.317");
    ++j;
  }
  if (resn == "ASN") {
    vector<int> NC(14);
    vector<int> BBC(14);
    vector<int> AT(14); /*ASN_ N  */
    AT[0] = 18;
    NC[0] = 2;
    BBC[0] = 0;
    /*ASN_ CA */ AT[1] = 19;
    NC[1] = 2;
    BBC[1] = 0;
    /*ASN_ C  */ AT[2] = 20;
    NC[2] = 2;
    BBC[2] = 0;
    /*ASN_ O  */ AT[3] = 21;
    NC[3] = 2;
    BBC[3] = 0;
    /*ASN_ CB */ AT[4] = 4;
    NC[4] = 2;
    BBC[4] = 2;
    /*ASN_ CG */ AT[5] = 1;
    NC[5] = 2;
    BBC[5] = 2;
    /*ASN_ OD1*/ AT[6] = 14;
    NC[6] = 2;
    BBC[6] = 2;
    /*ASN_ ND2*/ AT[7] = 9;
    NC[7] = 2;
    BBC[7] = 2;
    /*ASN_ H  */ AT[8] = 27;
    NC[8] = 2;
    BBC[8] = 0;
    /*ASN_ HA */ AT[9] = 25;
    NC[9] = 2;
    BBC[9] = 0;
    /*ASN_1HB */ AT[10] = 25;
    NC[10] = 2;
    BBC[10] = 1;
    /*ASN_2HB */ AT[11] = 25;
    NC[11] = 2;
    BBC[11] = 1;
    /*ASN_1HD2*/ AT[12] = 24;
    NC[12] = 2;
    BBC[12] = 2;
    /*ASN_2HD2*/ AT[13] = 24;
    NC[13] = 2;
    BBC[13] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   ASN     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  ASN     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   ASN     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   ASN     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  ASN     0       1.992  -0.780  -1.187");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  ASN     0       0.896  -1.345  -2.048");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  OD1 ASN     0      -0.293  -1.158  -1.763");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  ND2 ASN     0       1.272  -2.033  -3.095");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  H   ASN     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  HA  ASN     0       1.804  -0.481   0.916");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11 1HB  ASN     0       2.619  -0.129  -1.798");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 2HB  ASN     0       2.618  -1.599  -0.833");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 1HD2 ASN     0       0.586  -2.433  -3.704");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 2HD2 ASN     0       2.245  -2.159  -3.287");
    ++j;
  }
  if (resn == "PRO") {
    vector<int> NC(15);
    vector<int> BBC(15);
    vector<int> AT(15); /*PRO_ N  */
    AT[0] = 12;
    NC[0] = 3;
    BBC[0] = 0;
    /*PRO_ CA */ AT[1] = 19;
    NC[1] = 3;
    BBC[1] = 0;
    /*PRO_ C  */ AT[2] = 20;
    NC[2] = 3;
    BBC[2] = 2;
    /*PRO_ O  */ AT[3] = 21;
    NC[3] = 3;
    BBC[3] = 2;
    /*PRO_ CB */ AT[4] = 4;
    NC[4] = 3;
    BBC[4] = 2;
    /*PRO_ CG */ AT[5] = 4;
    NC[5] = 3;
    BBC[5] = 3;
    /*PRO_ CD */ AT[6] = 4;
    NC[6] = 3;
    BBC[6] = 3;
    /*PRO_ NV */ AT[7] = 57;
    NC[7] = 3;
    BBC[7] = 3;
    /*PRO_ HA */ AT[8] = 25;
    NC[8] = 3;
    BBC[8] = 2;
    /*PRO_1HB */ AT[9] = 25;
    NC[9] = 3;
    BBC[9] = 2;
    /*PRO_2HB */ AT[10] = 25;
    NC[10] = 3;
    BBC[10] = 2;
    /*PRO_1HG */ AT[11] = 25;
    NC[11] = 3;
    BBC[11] = 2;
    /*PRO_2HG */ AT[12] = 25;
    NC[12] = 3;
    BBC[12] = 2;
    /*PRO_1HD */ AT[13] = 25;
    NC[13] = 3;
    BBC[13] = 3;
    /*PRO_2HD */ AT[14] = 25;
    NC[14] = 3;
    BBC[14] = 3;
    int j = 0;
    if (true)
      *outiter++ =
          std::string("ATOM      1  N   PRO     0       0.000   0.000   0.000");
    ++j;
    if (true)
      *outiter++ =
          std::string("ATOM      2  CA  PRO     0       1.458   0.000   0.000");
    ++j;
    if (true)
      *outiter++ =
          std::string("ATOM      3  C   PRO     0       2.009   1.420   0.000");
    ++j;
    // if( true     ) *outiter++ = std::string("ATOM      4  O   PRO     0
    // 1.383   2.339  -0.529"); ++j;
    if (true)
      *outiter++ =
          std::string("ATOM      5  CB  PRO     0       1.803  -0.740  -1.297");
    ++j;
    if (true)
      *outiter++ =
          std::string("ATOM      6  CG  PRO     0       0.666  -0.428  -2.209");
    ++j;
    if (true)
      *outiter++ =
          std::string("ATOM      7  CD  PRO     0      -0.531  -0.340  -1.299");
    ++j;
    // if( true    ) *outiter++ = std::string("ATOM      8  NV  PRO     0
    // 0.075   0.007   0.001"); ++j;
    if (hydrogen)
      *outiter++ =
          std::string("ATOM      9  HA  PRO     0       1.820  -0.555   0.878");
    ++j;
    if (hydrogen)
      *outiter++ =
          std::string("ATOM     10 1HB  PRO     0       2.770  -0.388  -1.685");
    ++j;
    if (hydrogen)
      *outiter++ =
          std::string("ATOM     11 2HB  PRO     0       1.908  -1.816  -1.099");
    ++j;
    if (hydrogen)
      *outiter++ =
          std::string("ATOM     12 1HG  PRO     0       0.857   0.511  -2.748");
    ++j;
    if (hydrogen)
      *outiter++ =
          std::string("ATOM     13 2HG  PRO     0       0.562  -1.216  -2.970");
    ++j;
    if (hydrogen)
      *outiter++ =
          std::string("ATOM     14 1HD  PRO     0      -1.209   0.448  -1.659");
    ++j;
    if (hydrogen)
      *outiter++ =
          std::string("ATOM     15 2HD  PRO     0      -1.045  -1.312  -1.273");
    ++j;
  }
  if (resn == "GLN") {
    vector<int> NC(17);
    vector<int> BBC(17);
    vector<int> AT(17); /*GLN_ N  */
    AT[0] = 18;
    NC[0] = 3;
    BBC[0] = 0;
    /*GLN_ CA */ AT[1] = 19;
    NC[1] = 3;
    BBC[1] = 0;
    /*GLN_ C  */ AT[2] = 20;
    NC[2] = 3;
    BBC[2] = 0;
    /*GLN_ O  */ AT[3] = 21;
    NC[3] = 3;
    BBC[3] = 0;
    /*GLN_ CB */ AT[4] = 4;
    NC[4] = 3;
    BBC[4] = 2;
    /*GLN_ CG */ AT[5] = 4;
    NC[5] = 3;
    BBC[5] = 3;
    /*GLN_ CD */ AT[6] = 1;
    NC[6] = 3;
    BBC[6] = 3;
    /*GLN_ OE1*/ AT[7] = 14;
    NC[7] = 3;
    BBC[7] = 3;
    /*GLN_ NE2*/ AT[8] = 9;
    NC[8] = 3;
    BBC[8] = 3;
    /*GLN_ H  */ AT[9] = 27;
    NC[9] = 3;
    BBC[9] = 0;
    /*GLN_ HA */ AT[10] = 25;
    NC[10] = 3;
    BBC[10] = 0;
    /*GLN_1HB */ AT[11] = 25;
    NC[11] = 3;
    BBC[11] = 1;
    /*GLN_2HB */ AT[12] = 25;
    NC[12] = 3;
    BBC[12] = 1;
    /*GLN_1HG */ AT[13] = 25;
    NC[13] = 3;
    BBC[13] = 2;
    /*GLN_2HG */ AT[14] = 25;
    NC[14] = 3;
    BBC[14] = 2;
    /*GLN_1HE2*/ AT[15] = 24;
    NC[15] = 3;
    BBC[15] = 3;
    /*GLN_2HE2*/ AT[16] = 24;
    NC[16] = 3;
    BBC[16] = 3;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   GLN     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  GLN     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   GLN     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   GLN     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  GLN     0       1.994  -0.768  -1.211");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  GLN     0       0.914  -1.341  -2.113");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CD  GLN     0      -0.484  -1.025  -1.615");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  OE1 GLN     0      -0.657  -0.370  -0.583");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  NE2 GLN     0      -1.490  -1.489  -2.346");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  H   GLN     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  HA  GLN     0       1.804  -0.495   0.907");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 1HB  GLN     0       2.620  -0.109  -1.812");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 2HB  GLN     0       2.620  -1.593  -0.870");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 1HG  GLN     0       1.025  -0.916  -3.111");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 2HG  GLN     0       1.025  -2.424  -2.153");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16 1HE2 GLN     0      -2.435  -1.312  -2.067");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17 2HE2 GLN     0      -1.305  -2.015  -3.176");
    ++j;
  }
  if (resn == "ARG") {
    vector<int> NC(24);
    vector<int> BBC(24);
    vector<int> AT(24); /*ARG_ N  */
    AT[0] = 18;
    NC[0] = 4;
    BBC[0] = 0;
    /*ARG_ CA */ AT[1] = 19;
    NC[1] = 4;
    BBC[1] = 0;
    /*ARG_ C  */ AT[2] = 20;
    NC[2] = 4;
    BBC[2] = 2;
    /*ARG_ O  */ AT[3] = 21;
    NC[3] = 4;
    BBC[3] = 2;
    /*ARG_ CB */ AT[4] = 4;
    NC[4] = 4;
    BBC[4] = 2;
    /*ARG_ CG */ AT[5] = 4;
    NC[5] = 4;
    BBC[5] = 3;
    /*ARG_ CD */ AT[6] = 4;
    NC[6] = 4;
    BBC[6] = 4;
    /*ARG_ NE */ AT[7] = 11;
    NC[7] = 4;
    BBC[7] = 4;
    /*ARG_ CZ */ AT[8] = 6;
    NC[8] = 4;
    BBC[8] = 4;
    /*ARG_ NH1*/ AT[9] = 11;
    NC[9] = 4;
    BBC[9] = 4;
    /*ARG_ NH2*/ AT[10] = 11;
    NC[10] = 4;
    BBC[10] = 4;
    /*ARG_ H  */ AT[11] = 27;
    NC[11] = 4;
    BBC[11] = 2;
    /*ARG_ HA */ AT[12] = 25;
    NC[12] = 4;
    BBC[12] = 2;
    /*ARG_1HB */ AT[13] = 25;
    NC[13] = 4;
    BBC[13] = 2;
    /*ARG_2HB */ AT[14] = 25;
    NC[14] = 4;
    BBC[14] = 2;
    /*ARG_1HG */ AT[15] = 25;
    NC[15] = 4;
    BBC[15] = 2;
    /*ARG_2HG */ AT[16] = 25;
    NC[16] = 4;
    BBC[16] = 2;
    /*ARG_1HD */ AT[17] = 25;
    NC[17] = 4;
    BBC[17] = 3;
    /*ARG_2HD */ AT[18] = 25;
    NC[18] = 4;
    BBC[18] = 3;
    /*ARG_ HE */ AT[19] = 24;
    NC[19] = 4;
    BBC[19] = 4;
    /*ARG_1HH1*/ AT[20] = 24;
    NC[20] = 4;
    BBC[20] = 4;
    /*ARG_2HH1*/ AT[21] = 24;
    NC[21] = 4;
    BBC[21] = 4;
    /*ARG_1HH2*/ AT[22] = 24;
    NC[22] = 4;
    BBC[22] = 4;
    /*ARG_2HH2*/ AT[23] = 24;
    NC[23] = 4;
    BBC[23] = 4;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   ARG     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  ARG     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   ARG     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   ARG     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  ARG     0       1.993  -0.748  -1.212");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  ARG     0       0.927  -1.318  -2.134");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CD  ARG     0      -0.442  -1.015  -1.644");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  NE  ARG     0      -0.418  -0.251  -0.407");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  CZ  ARG     0       0.703   0.141   0.229");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  NH1 ARG     0       1.881  -0.165  -0.267");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  NH2 ARG     0       0.618   0.834   1.351");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12  H   ARG     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13  HA  ARG     0       1.804  -0.506   0.902");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 1HB  ARG     0       2.619  -0.082  -1.804");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 2HB  ARG     0       2.620  -1.577  -0.880");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16 1HG  ARG     0       1.038  -0.886  -3.129");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17 2HG  ARG     0       1.038  -2.401  -2.194");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     18 1HD  ARG     0      -0.976  -0.433  -2.394");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     19 2HD  ARG     0      -0.976  -1.946  -1.459");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     20  HE  ARG     0      -1.306   0.003   0.005");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     21 1HH1 ARG     0       1.946  -0.695  -1.125");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     22 2HH1 ARG     0       2.721   0.129   0.209");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     23 1HH2 ARG     0      -0.288   1.070   1.732");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     24 2HH2 ARG     0       1.458   1.128   1.827");
    ++j;
  }
  if (resn == "SER") {
    vector<int> NC(11);
    vector<int> BBC(11);
    vector<int> AT(11); /*SER_ N  */
    AT[0] = 18;
    NC[0] = 1;
    BBC[0] = 0;
    /*SER_ CA */ AT[1] = 19;
    NC[1] = 1;
    BBC[1] = 2;
    /*SER_ C  */ AT[2] = 20;
    NC[2] = 1;
    BBC[2] = 0;
    /*SER_ O  */ AT[3] = 21;
    NC[3] = 1;
    BBC[3] = 0;
    /*SER_ CB */ AT[4] = 4;
    NC[4] = 1;
    BBC[4] = 2;
    /*SER_ OG */ AT[5] = 13;
    NC[5] = 1;
    BBC[5] = 2;
    /*SER_ H  */ AT[6] = 27;
    NC[6] = 1;
    BBC[6] = 0;
    /*SER_ HA */ AT[7] = 25;
    NC[7] = 1;
    BBC[7] = 0;
    /*SER_1HB */ AT[8] = 25;
    NC[8] = 1;
    BBC[8] = 1;
    /*SER_2HB */ AT[9] = 25;
    NC[9] = 1;
    BBC[9] = 1;
    /*SER_ HG */ AT[10] = 24;
    NC[10] = 1;
    BBC[10] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   SER     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  SER     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   SER     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   SER     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  SER     0       1.980  -0.754  -1.207");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  OG  SER     0       0.925  -1.243  -1.989");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  H   SER     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  HA  SER     0       1.804  -0.502   0.905");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9 1HB  SER     0       2.605  -0.093  -1.807");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10 2HB  SER     0       2.605  -1.583  -0.876");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  HG  SER     0       0.118  -0.967  -1.547");
    ++j;
  }
  if (resn == "THR") {
    vector<int> NC(14);
    vector<int> BBC(14);
    vector<int> AT(14); /*THR_ N  */
    AT[0] = 18;
    NC[0] = 1;
    BBC[0] = 0;
    /*THR_ CA */ AT[1] = 19;
    NC[1] = 1;
    BBC[1] = 2;
    /*THR_ C  */ AT[2] = 20;
    NC[2] = 1;
    BBC[2] = 0;
    /*THR_ O  */ AT[3] = 21;
    NC[3] = 1;
    BBC[3] = 0;
    /*THR_ CB */ AT[4] = 3;
    NC[4] = 1;
    BBC[4] = 2;
    /*THR_ OG1*/ AT[5] = 13;
    NC[5] = 1;
    BBC[5] = 2;
    /*THR_ CG2*/ AT[6] = 5;
    NC[6] = 1;
    BBC[6] = 1;
    /*THR_ H  */ AT[7] = 27;
    NC[7] = 1;
    BBC[7] = 0;
    /*THR_ HA */ AT[8] = 25;
    NC[8] = 1;
    BBC[8] = 0;
    /*THR_ HB */ AT[9] = 25;
    NC[9] = 1;
    BBC[9] = 1;
    /*THR_ HG1*/ AT[10] = 24;
    NC[10] = 1;
    BBC[10] = 2;
    /*THR_1HG2*/ AT[11] = 25;
    NC[11] = 1;
    BBC[11] = 1;
    /*THR_2HG2*/ AT[12] = 25;
    NC[12] = 1;
    BBC[12] = 1;
    /*THR_3HG2*/ AT[13] = 25;
    NC[13] = 1;
    BBC[13] = 1;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   THR     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  THR     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   THR     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   THR     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  THR     0       2.012  -0.770  -1.213");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  OG1 THR     0       0.925  -1.271  -2.002");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CG2 THR     0       2.880  -1.933  -0.756");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  H   THR     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  HA  THR     0       1.804  -0.495   0.908");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  HB  THR     0       2.610  -0.098  -1.828");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  HG1 THR     0       0.093  -1.014  -1.598");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 1HG2 THR     0       3.263  -2.465  -1.626");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 2HG2 THR     0       3.714  -1.554  -0.165");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 3HG2 THR     0       2.285  -2.613  -0.148");
    ++j;
  }
  if (resn == "VAL") {
    vector<int> NC(16);
    vector<int> BBC(16);
    vector<int> AT(16); /*VAL_ N  */
    AT[0] = 18;
    NC[0] = 1;
    BBC[0] = 0;
    /*VAL_ CA */ AT[1] = 19;
    NC[1] = 1;
    BBC[1] = 1;
    /*VAL_ C  */ AT[2] = 20;
    NC[2] = 1;
    BBC[2] = 0;
    /*VAL_ O  */ AT[3] = 21;
    NC[3] = 1;
    BBC[3] = 0;
    /*VAL_ CB */ AT[4] = 3;
    NC[4] = 1;
    BBC[4] = 1;
    /*VAL_ CG1*/ AT[5] = 5;
    NC[5] = 1;
    BBC[5] = 1;
    /*VAL_ CG2*/ AT[6] = 5;
    NC[6] = 1;
    BBC[6] = 1;
    /*VAL_ H  */ AT[7] = 27;
    NC[7] = 1;
    BBC[7] = 0;
    /*VAL_ HA */ AT[8] = 25;
    NC[8] = 1;
    BBC[8] = 0;
    /*VAL_ HB */ AT[9] = 25;
    NC[9] = 1;
    BBC[9] = 1;
    /*VAL_1HG1*/ AT[10] = 25;
    NC[10] = 1;
    BBC[10] = 1;
    /*VAL_2HG1*/ AT[11] = 25;
    NC[11] = 1;
    BBC[11] = 1;
    /*VAL_3HG1*/ AT[12] = 25;
    NC[12] = 1;
    BBC[12] = 1;
    /*VAL_1HG2*/ AT[13] = 25;
    NC[13] = 1;
    BBC[13] = 1;
    /*VAL_2HG2*/ AT[14] = 25;
    NC[14] = 1;
    BBC[14] = 1;
    /*VAL_3HG2*/ AT[15] = 25;
    NC[15] = 1;
    BBC[15] = 1;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   VAL     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  VAL     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   VAL     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   VAL     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  VAL     0       1.992  -0.755  -1.232");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG1 VAL     0       0.840  -1.274  -2.079");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CG2 VAL     0       2.894   0.159  -2.048");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  H   VAL     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  HA  VAL     0       1.804  -0.509   0.900");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  HB  VAL     0       2.561  -1.622  -0.896");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11 1HG1 VAL     0       1.236  -1.805  -2.945");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12 2HG1 VAL     0       0.230  -1.954  -1.485");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13 3HG1 VAL     0       0.230  -0.436  -2.416");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14 1HG2 VAL     0       3.268  -0.381  -2.917");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 2HG2 VAL     0       2.326   1.029  -2.379");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16 3HG2 VAL     0       3.733   0.485  -1.434");
    ++j;
  }
  if (resn == "TRP") {
    vector<int> NC(24);
    vector<int> BBC(24);
    vector<int> AT(24); /*TRP_ N  */
    AT[0] = 18;
    NC[0] = 2;
    BBC[0] = 0;
    /*TRP_ CA */ AT[1] = 19;
    NC[1] = 2;
    BBC[1] = 0;
    /*TRP_ C  */ AT[2] = 20;
    NC[2] = 2;
    BBC[2] = 0;
    /*TRP_ O  */ AT[3] = 21;
    NC[3] = 2;
    BBC[3] = 0;
    /*TRP_ CB */ AT[4] = 4;
    NC[4] = 2;
    BBC[4] = 2;
    /*TRP_ CG */ AT[5] = 6;
    NC[5] = 2;
    BBC[5] = 2;
    /*TRP_ CD1*/ AT[6] = 6;
    NC[6] = 2;
    BBC[6] = 2;
    /*TRP_ CD2*/ AT[7] = 6;
    NC[7] = 2;
    BBC[7] = 2;
    /*TRP_ NE1*/ AT[8] = 7;
    NC[8] = 2;
    BBC[8] = 2;
    /*TRP_ CE2*/ AT[9] = 6;
    NC[9] = 2;
    BBC[9] = 2;
    /*TRP_ CE3*/ AT[10] = 6;
    NC[10] = 2;
    BBC[10] = 2;
    /*TRP_ CZ2*/ AT[11] = 6;
    NC[11] = 2;
    BBC[11] = 2;
    /*TRP_ CZ3*/ AT[12] = 6;
    NC[12] = 2;
    BBC[12] = 2;
    /*TRP_ CH2*/ AT[13] = 6;
    NC[13] = 2;
    BBC[13] = 2;
    /*TRP_ H  */ AT[14] = 27;
    NC[14] = 2;
    BBC[14] = 0;
    /*TRP_ HA */ AT[15] = 25;
    NC[15] = 2;
    BBC[15] = 0;
    /*TRP_1HB */ AT[16] = 25;
    NC[16] = 2;
    BBC[16] = 1;
    /*TRP_2HB */ AT[17] = 25;
    NC[17] = 2;
    BBC[17] = 1;
    /*TRP_ HD1*/ AT[18] = 26;
    NC[18] = 2;
    BBC[18] = 2;
    /*TRP_ HE1*/ AT[19] = 24;
    NC[19] = 2;
    BBC[19] = 2;
    /*TRP_ HE3*/ AT[20] = 26;
    NC[20] = 2;
    BBC[20] = 2;
    /*TRP_ HZ2*/ AT[21] = 26;
    NC[21] = 2;
    BBC[21] = 2;
    /*TRP_ HZ3*/ AT[22] = 26;
    NC[22] = 2;
    BBC[22] = 2;
    /*TRP_ HH2*/ AT[23] = 26;
    NC[23] = 2;
    BBC[23] = 2;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   TRP     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  TRP     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   TRP     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   TRP     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  TRP     0       1.991  -0.758  -1.217");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  TRP     0       0.912  -1.307  -2.100");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CD1 TRP     0      -0.433  -1.190  -1.913");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  CD2 TRP     0       1.082  -2.067  -3.322");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  NE1 TRP     0      -1.109  -1.823  -2.927");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  CE2 TRP     0      -0.197  -2.363  -3.799");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  CE3 TRP     0       2.201  -2.511  -4.037");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12  CZ2 TRP     0      -0.396  -3.088  -4.963");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13  CZ3 TRP     0       2.001  -3.237  -5.205");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14  CH2 TRP     0       0.736  -3.517  -5.656");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15  H   TRP     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16  HA  TRP     0       1.804  -0.503   0.903");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17 1HB  TRP     0       2.617  -0.094  -1.814");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     18 2HB  TRP     0       2.617  -1.586  -0.885");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     19  HD1 TRP     0      -0.903  -0.671  -1.079");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     20  HE1 TRP     0      -2.113  -1.882  -3.017");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     21  HE3 TRP     0       3.208  -2.289  -3.685");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     22  HZ2 TRP     0      -1.393  -3.321  -5.337");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     23  HZ3 TRP     0       2.878  -3.580  -5.756");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     24  HH2 TRP     0       0.616  -4.089  -6.577");
    ++j;
  }
  if (resn == "TYR") {
    vector<int> NC(21);
    vector<int> BBC(21);
    vector<int> AT(21); /*TYR_ N  */
    AT[0] = 18;
    NC[0] = 2;
    BBC[0] = 0;
    /*TYR_ CA */ AT[1] = 19;
    NC[1] = 2;
    BBC[1] = 0;
    /*TYR_ C  */ AT[2] = 20;
    NC[2] = 2;
    BBC[2] = 0;
    /*TYR_ O  */ AT[3] = 21;
    NC[3] = 2;
    BBC[3] = 0;
    /*TYR_ CB */ AT[4] = 4;
    NC[4] = 2;
    BBC[4] = 2;
    /*TYR_ CG */ AT[5] = 6;
    NC[5] = 2;
    BBC[5] = 2;
    /*TYR_ CD1*/ AT[6] = 6;
    NC[6] = 2;
    BBC[6] = 2;
    /*TYR_ CD2*/ AT[7] = 6;
    NC[7] = 2;
    BBC[7] = 2;
    /*TYR_ CE1*/ AT[8] = 6;
    NC[8] = 2;
    BBC[8] = 2;
    /*TYR_ CE2*/ AT[9] = 6;
    NC[9] = 2;
    BBC[9] = 2;
    /*TYR_ CZ */ AT[10] = 6;
    NC[10] = 2;
    BBC[10] = 3;
    /*TYR_ OH */ AT[11] = 13;
    NC[11] = 2;
    BBC[11] = 3;
    /*TYR_ H  */ AT[12] = 27;
    NC[12] = 2;
    BBC[12] = 0;
    /*TYR_ HA */ AT[13] = 25;
    NC[13] = 2;
    BBC[13] = 0;
    /*TYR_1HB */ AT[14] = 25;
    NC[14] = 2;
    BBC[14] = 1;
    /*TYR_2HB */ AT[15] = 25;
    NC[15] = 2;
    BBC[15] = 1;
    /*TYR_ HD1*/ AT[16] = 26;
    NC[16] = 2;
    BBC[16] = 2;
    /*TYR_ HD2*/ AT[17] = 26;
    NC[17] = 2;
    BBC[17] = 2;
    /*TYR_ HE1*/ AT[18] = 26;
    NC[18] = 2;
    BBC[18] = 2;
    /*TYR_ HE2*/ AT[19] = 26;
    NC[19] = 2;
    BBC[19] = 2;
    /*TYR_ HH */ AT[20] = 24;
    NC[20] = 2;
    BBC[20] = 3;
    int j = 0;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      1  N   TYR     0       0.000   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      2  CA  TYR     0       1.458   0.000   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      3  C   TYR     0       2.009   1.420   0.000");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      4  O   TYR     0       1.383   2.339  -0.529");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      5  CB  TYR     0       1.993  -0.777  -1.205");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      6  CG  TYR     0       0.910  -1.349  -2.093");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      7  CD1 TYR     0      -0.425  -1.145  -1.777");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      8  CD2 TYR     0       1.252  -2.077  -3.222");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM      9  CE1 TYR     0      -1.414  -1.667  -2.587");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     10  CE2 TYR     0       0.263  -2.599  -4.033");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     11  CZ  TYR     0      -1.065  -2.397  -3.719");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     12  OH  TYR     0      -2.050  -2.917  -4.526");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     13  H   TYR     0      -0.492   0.764  -0.441");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     14  HA  TYR     0       1.804  -0.489   0.911");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     15 1HB  TYR     0       2.620  -0.121  -1.812");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     16 2HB  TYR     0       2.619  -1.599  -0.859");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     17  HD1 TYR     0      -0.694  -0.573  -0.889");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     18  HD2 TYR     0       2.301  -2.238  -3.470");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     19  HE1 TYR     0      -2.463  -1.507  -2.339");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     20  HE2 TYR     0       0.532  -3.172  -4.921");
    ++j;
    if (check(AT[j], NC[j], BBC[j]))
      *outiter++ =
          std::string("ATOM     21  HH  TYR     0      -1.646  -3.389  -5.258");
    ++j;
  }

  return lines;
}
}
}  // end namespace
