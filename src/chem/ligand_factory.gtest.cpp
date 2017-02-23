#include <gtest/gtest.h>

#include <fstream>
#include <iterator>  // std::back_inserter
#include "actor/Atom.hpp"
#include "ligand_factory.hpp"
#include "util/SimpleArray.hpp"

#include <boost/foreach.hpp>

namespace scheme {
namespace chemical {
namespace test {

using std::cout;
using std::endl;

TEST(ligand_factory, make_atom_pdbline) {
  typedef util::SimpleArray<3, double> Position;
  typedef actor::Atom<Position> Atom;
  LigandFactory<Atom> f;

  std::string l;
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(false, f.make_atom_pdbline(l).data().ishet);
  l = "HETATM    7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(true, f.make_atom_pdbline(l).data().ishet);

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(7, f.make_atom_pdbline(l).data().atomnum);
  l = "ATOM  999999 C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(999999, f.make_atom_pdbline(l).data().atomnum);
  l = "ATOM  9999999C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(999999, f.make_atom_pdbline(l).data().atomnum);

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ("C3", f.make_atom_pdbline(l).data().atomname);
  l = "ATOM      7 ATOMABTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ("ATOM", f.make_atom_pdbline(l).data().atomname);

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ("BTN", f.make_atom_pdbline(l).data().resname);
  l = "ATOM      7 ATOMABTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ("ABTN", f.make_atom_pdbline(l).data().resname);

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ('X', f.make_atom_pdbline(l).data().chain);
  l = "ATOM      7  C3  BTN     1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(' ', f.make_atom_pdbline(l).data().chain);

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(1, f.make_atom_pdbline(l).data().resnum);
  l = "ATOM      7  C3  BTN X9999      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(9999, f.make_atom_pdbline(l).data().resnum);
  l = "ATOM      7  C3  BTN X9999      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(9999, f.make_atom_pdbline(l).data().resnum);
  l = "ATOM      7  C3  BTN X999999    -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(999999, f.make_atom_pdbline(l).data().resnum);

  l = "ATOM      7  C3  BTN X   1  -99999.999 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(-99999.999, f.make_atom_pdbline(l).position()[0]);
  l = "ATOM      7  C3  BTN X   1      -0.470-999.999   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(-999.999, f.make_atom_pdbline(l).position()[1]);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087-999.999  1.00 20.00      "
      "     C";
  ASSERT_EQ(-999.999, f.make_atom_pdbline(l).position()[2]);

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377111.00 20.00      "
      "     C";
  ASSERT_EQ(111, f.make_atom_pdbline(l).data().occ);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00120.00      "
      "     C";
  ASSERT_EQ(120, f.make_atom_pdbline(l).data().bfac);

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ("C", f.make_atom_pdbline(l).data().elem);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00  "
      "AAAAAAAAAA";
  ASSERT_EQ("AAAAAAAAAA", f.make_atom_pdbline(l).data().elem);

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(l, io::dump_pdb_atom(f.make_atom_pdbline(l)));
  l = "ATOM      7  C   BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(l, io::dump_pdb_atom(f.make_atom_pdbline(l)));
  l = "ATOM      7  C3A BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(l, io::dump_pdb_atom(f.make_atom_pdbline(l)));
  l = "ATOM      7 AC3A BTN X   1      -0.470 -12.087   5.377  1.00 20.00      "
      "     C";
  ASSERT_EQ(l, io::dump_pdb_atom(f.make_atom_pdbline(l)));

  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00";
  ASSERT_EQ(false, f.make_atom_pdbline(l).data().ishet);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00";
  ASSERT_EQ(7, f.make_atom_pdbline(l).data().atomnum);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00";
  ASSERT_EQ("C3", f.make_atom_pdbline(l).data().atomname);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00";
  ASSERT_EQ("BTN", f.make_atom_pdbline(l).data().resname);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00";
  ASSERT_EQ('X', f.make_atom_pdbline(l).data().chain);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00";
  ASSERT_EQ(1, f.make_atom_pdbline(l).data().resnum);
  l = "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 20.00";
  ASSERT_EQ("", f.make_atom_pdbline(l).data().elem);
}

TEST(ligand_factory, make_biotin) {
  typedef util::SimpleArray<3, double> Position;
  typedef actor::Atom<Position> Atom;
  LigandFactory<Atom> f;

  std::vector<Atom> btn;
  f.make_biotin_minimal(std::back_inserter(btn));
  // BOOST_FOREACH(Atom a,btn) cout << io::dump_pdb_atom(a) << endl;
  ASSERT_EQ(io::dump_pdb_atom(btn[0]),
            "ATOM      1  N1  BTN X   1       0.696 -12.422   3.375  1.00 "
            "20.00           N");
  ASSERT_EQ(io::dump_pdb_atom(btn[1]),
            "ATOM      2  S1  BTN X   1       0.576  -9.666   5.336  1.00 "
            "20.00           S");
  ASSERT_EQ(io::dump_pdb_atom(btn[2]),
            "ATOM      3  C1  BTN X   1      -0.523 -10.824   6.189  1.00 "
            "20.00           C");
  ASSERT_EQ(io::dump_pdb_atom(btn[3]),
            "ATOM      4  N2  BTN X   1      -1.324 -12.123   4.201  1.00 "
            "20.00           N");
  ASSERT_EQ(io::dump_pdb_atom(btn[4]),
            "ATOM      5  C2  BTN X   1      -0.608 -12.327   3.072  1.00 "
            "20.00           C");
  ASSERT_EQ(io::dump_pdb_atom(btn[5]),
            "ATOM      6  O1  BTN X   1      -1.125 -12.422   1.933  1.00 "
            "20.00           O");
  ASSERT_EQ(io::dump_pdb_atom(btn[6]),
            "ATOM      7  C3  BTN X   1      -0.470 -12.087   5.377  1.00 "
            "20.00           C");
  ASSERT_EQ(io::dump_pdb_atom(btn[7]),
            "ATOM      8  C4  BTN X   1       0.953 -12.267   4.780  1.00 "
            "20.00           C");
  ASSERT_EQ(io::dump_pdb_atom(btn[8]),
            "ATOM      9  C5  BTN X   1       1.765 -11.040   5.134  1.00 "
            "20.00           C");
  ASSERT_EQ(io::dump_pdb_atom(btn[9]),
            "ATOM     10  C6  BTN X   1      -1.836 -10.395   6.850  1.00 "
            "20.00           C");
}

TEST(ligand_factory, make_aas) {
  typedef util::SimpleArray<3, double> Position;
  typedef actor::Atom<Position> Atom;
  LigandFactory<Atom> f;

  std::vector<Atom> GLY;
  f.make_atoms(std::back_inserter(GLY), "GLY", false);
  EXPECT_EQ(GLY.size(), 3);
  BOOST_FOREACH (Atom a, GLY)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,GLY)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> GLY_H;
  f.make_atoms(std::back_inserter(GLY_H), "GLY", true);
  EXPECT_EQ(GLY_H.size(), 5);
  BOOST_FOREACH (Atom a, GLY_H)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,GLY_H)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> ALA;
  f.make_atoms(std::back_inserter(ALA), "ALA", false);
  EXPECT_EQ(ALA.size(), 4);
  BOOST_FOREACH (Atom a, ALA)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ALA)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> ALA_H;
  f.make_atoms(std::back_inserter(ALA_H), "ALA", true);
  EXPECT_EQ(ALA_H.size(), 8);
  BOOST_FOREACH (Atom a, ALA_H)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ALA_H)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> CYS;
  f.make_atoms(std::back_inserter(CYS), "CYS", false);
  EXPECT_EQ(CYS.size(), 3);
  BOOST_FOREACH (Atom a, CYS)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,CYS)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> CYS_H;
  f.make_atoms(std::back_inserter(CYS_H), "CYS", true);
  EXPECT_EQ(CYS_H.size(), 5);
  BOOST_FOREACH (Atom a, CYS_H)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,CYS_H)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> PHE;
  f.make_atoms(std::back_inserter(PHE), "PHE", false);
  EXPECT_EQ(PHE.size(), 7);
  BOOST_FOREACH (Atom a, PHE)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,PHE)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> PHE_H;
  f.make_atoms(std::back_inserter(PHE_H), "PHE", true);
  EXPECT_EQ(PHE_H.size(), 12);
  BOOST_FOREACH (Atom a, PHE_H)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,PHE_H)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> TYR;
  f.make_atoms(std::back_inserter(TYR), "TYR", false);
  EXPECT_EQ(TYR.size(), 8);
  BOOST_FOREACH (Atom a, TYR)
    ASSERT_GT(a.type(), 0);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,TYR)
  // io::dump_pdb_atom(out,a); out.close();

  std::vector<Atom> TYR_H;
  f.make_atoms(std::back_inserter(TYR_H), "TYR", true);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,TYR_H)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,TYR_H) cout << a << endl;
  EXPECT_EQ(TYR_H.size(), 12);
  BOOST_FOREACH (Atom a, TYR_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> ASP_H;
  f.make_atoms(std::back_inserter(ASP_H), "ASP", true);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ASP_H)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,ASP_H) cout << a << endl;
  EXPECT_EQ(ASP_H.size(), 4);
  BOOST_FOREACH (Atom a, ASP_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> GLU_H;
  f.make_atoms(std::back_inserter(GLU_H), "GLU", true);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,GLU_H)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,GLU_H) cout << a << endl;
  EXPECT_EQ(GLU_H.size(), 4);
  BOOST_FOREACH (Atom a, GLU_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> ASN_H;
  f.make_atoms(std::back_inserter(ASN_H), "ASN", true);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ASN_H)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,ASN_H) cout << a << endl;
  EXPECT_EQ(ASN_H.size(), 6);
  BOOST_FOREACH (Atom a, ASN_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> ASN;
  f.make_atoms(std::back_inserter(ASN), "ASN", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ASN)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,ASN) cout << a << endl;
  EXPECT_EQ(ASN.size(), 4);
  BOOST_FOREACH (Atom a, ASN)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> GLN_H;
  f.make_atoms(std::back_inserter(GLN_H), "GLN", true);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,GLN_H)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,GLN_H) cout << a << endl;
  EXPECT_EQ(GLN_H.size(), 6);
  BOOST_FOREACH (Atom a, GLN_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> GLN;
  f.make_atoms(std::back_inserter(GLN), "GLN", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,GLN)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,GLN) cout << a << endl;
  EXPECT_EQ(GLN.size(), 4);
  BOOST_FOREACH (Atom a, GLN)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> TRP_H;
  f.make_atoms(std::back_inserter(TRP_H), "TRP", true);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,TRP_H)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,TRP_H) cout << a << endl;
  EXPECT_EQ(TRP_H.size(), 16);
  BOOST_FOREACH (Atom a, TRP_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> TRP;
  f.make_atoms(std::back_inserter(TRP), "TRP", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,TRP)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,TRP) cout << a << endl;
  EXPECT_EQ(TRP.size(), 10);
  BOOST_FOREACH (Atom a, TRP)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> LEU_H;
  f.make_atoms(std::back_inserter(LEU_H), "LEU", true);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,LEU_H)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,LEU_H) cout << a << endl;
  EXPECT_EQ(LEU_H.size(), 11);
  BOOST_FOREACH (Atom a, LEU_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> LEU;
  f.make_atoms(std::back_inserter(LEU), "LEU", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,LEU)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,LEU) cout << a << endl;
  EXPECT_EQ(LEU.size(), 4);
  BOOST_FOREACH (Atom a, LEU)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> VAL_H;
  f.make_atoms(std::back_inserter(VAL_H), "VAL", true);
  // std::ofstream vhout("testval.pdb"); BOOST_FOREACH(Atom a,VAL_H)
  // io::dump_pdb_atom(vhout,a); vhout.close();
  // BOOST_FOREACH(Atom a,VAL_H) cout << a << endl;
  EXPECT_EQ(VAL_H.size(), 11);
  BOOST_FOREACH (Atom a, VAL_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> VAL;
  f.make_atoms(std::back_inserter(VAL), "VAL", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,VAL)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,VAL) cout << a << endl;
  EXPECT_EQ(VAL.size(), 4);
  BOOST_FOREACH (Atom a, VAL)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> MET_H;
  f.make_atoms(std::back_inserter(MET_H), "MET", true);
  // std::ofstream mhout("test.pdb"); BOOST_FOREACH(Atom a,MET_H)
  // io::dump_pdb_atom(mhout,a); mhout.close();
  // BOOST_FOREACH(Atom a,MET_H) cout << a << endl;
  EXPECT_EQ(MET_H.size(), 6);
  BOOST_FOREACH (Atom a, MET_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> MET;
  f.make_atoms(std::back_inserter(MET), "MET", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,MET)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,MET) cout << a << endl;
  EXPECT_EQ(MET.size(), 3);
  BOOST_FOREACH (Atom a, MET)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> PRO_H;
  f.make_atoms(std::back_inserter(PRO_H), "PRO", true);
  // std::ofstream mhout("test.pdb"); BOOST_FOREACH(Atom a,PRO_H)
  // io::dump_pdb_atom(mhout,a); mhout.close();
  // BOOST_FOREACH(Atom a,PRO_H) cout << a << endl;
  EXPECT_EQ(PRO_H.size(), 13);
  BOOST_FOREACH (Atom a, PRO_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> PRO;
  f.make_atoms(std::back_inserter(PRO), "PRO", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,PRO)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,PRO) cout << a << endl;
  EXPECT_EQ(PRO.size(), 6);
  BOOST_FOREACH (Atom a, PRO)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> HIS_H;
  f.make_atoms(std::back_inserter(HIS_H), "HIS", true);
  // std::ofstream mhout("test.pdb"); BOOST_FOREACH(Atom a,HIS_H)
  // io::dump_pdb_atom(mhout,a); mhout.close();
  // BOOST_FOREACH(Atom a,HIS_H) cout << a << endl;
  EXPECT_EQ(HIS_H.size(), 8);
  BOOST_FOREACH (Atom a, HIS_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> HIS;
  f.make_atoms(std::back_inserter(HIS), "HIS", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,HIS)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,HIS) cout << a << endl;
  EXPECT_EQ(HIS.size(), 6);
  BOOST_FOREACH (Atom a, HIS)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> ARG_H;
  f.make_atoms(std::back_inserter(ARG_H), "ARG", true);
  // std::ofstream mhout("test.pdb"); BOOST_FOREACH(Atom a,ARG_H)
  // io::dump_pdb_atom(mhout,a); mhout.close();
  // BOOST_FOREACH(Atom a,ARG_H) cout << a << endl;
  EXPECT_EQ(ARG_H.size(), 10);
  BOOST_FOREACH (Atom a, ARG_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> ARG;
  f.make_atoms(std::back_inserter(ARG), "ARG", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ARG)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,ARG) cout << a << endl;
  EXPECT_EQ(ARG.size(), 5);
  BOOST_FOREACH (Atom a, ARG)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> LYS_H;
  f.make_atoms(std::back_inserter(LYS_H), "LYS", true);
  // std::ofstream mhout("test.pdb"); BOOST_FOREACH(Atom a,LYS_H)
  // io::dump_pdb_atom(mhout,a); mhout.close();
  // BOOST_FOREACH(Atom a,LYS_H) cout << a << endl;
  EXPECT_EQ(LYS_H.size(), 8);
  BOOST_FOREACH (Atom a, LYS_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> LYS;
  f.make_atoms(std::back_inserter(LYS), "LYS", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,LYS)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,LYS) cout << a << endl;
  EXPECT_EQ(LYS.size(), 3);
  BOOST_FOREACH (Atom a, LYS)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> ILE_H;
  f.make_atoms(std::back_inserter(ILE_H), "ILE", true);
  // std::ofstream mhout("test.pdb"); BOOST_FOREACH(Atom a,ILE_H)
  // io::dump_pdb_atom(mhout,a); mhout.close();
  // BOOST_FOREACH(Atom a,ILE_H) cout << a << endl;
  EXPECT_EQ(ILE_H.size(), 8);
  BOOST_FOREACH (Atom a, ILE_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> ILE;
  f.make_atoms(std::back_inserter(ILE), "ILE", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,ILE)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,ILE) cout << a << endl;
  EXPECT_EQ(ILE.size(), 3);
  BOOST_FOREACH (Atom a, ILE)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> SER_H;
  f.make_atoms(std::back_inserter(SER_H), "SER", true);
  // std::ofstream shout("test.pdb"); BOOST_FOREACH(Atom a,SER_H)
  // io::dump_pdb_atom(shout,a); shout.close();
  // BOOST_FOREACH(Atom a,SER_H) cout << a << endl;
  EXPECT_EQ(SER_H.size(), 5);
  BOOST_FOREACH (Atom a, SER_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> SER;
  f.make_atoms(std::back_inserter(SER), "SER", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,SER)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,SER) cout << a << endl;
  EXPECT_EQ(SER.size(), 3);
  BOOST_FOREACH (Atom a, SER)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> THR_H;
  f.make_atoms(std::back_inserter(THR_H), "THR", true);
  // std::ofstream shout("test.pdb"); BOOST_FOREACH(Atom a,THR_H)
  // io::dump_pdb_atom(shout,a); shout.close();
  // BOOST_FOREACH(Atom a,THR_H) cout << a << endl;
  EXPECT_EQ(THR_H.size(), 8);
  BOOST_FOREACH (Atom a, THR_H)
    ASSERT_GT(a.type(), 0);

  std::vector<Atom> THR;
  f.make_atoms(std::back_inserter(THR), "THR", false);
  // std::ofstream out("test.pdb"); BOOST_FOREACH(Atom a,THR)
  // io::dump_pdb_atom(out,a); out.close();
  // BOOST_FOREACH(Atom a,THR) cout << a << endl;
  EXPECT_EQ(THR.size(), 4);
  BOOST_FOREACH (Atom a, THR)
    ASSERT_GT(a.type(), 0);
}
}
}
}
