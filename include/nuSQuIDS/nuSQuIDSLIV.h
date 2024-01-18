#ifndef nusquidslv_H
#define nusquidslv_H

#include <vector>
#include <iostream>
#include <nuSQuIDS/nuSQuIDS.h>
#include <iomanip>
#include <math.h>
#include <cmath>


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>



// // Printing matrix function
void print_gsl_matrix(gsl_matrix_complex* matrix) {
  for (size_t i = 0; i < matrix->size1; i++) {
    for (size_t j = 0; j < matrix->size2; j++) {
      //std::cout << std::fixed << std::setprecision(2) << GSL_REAL(gsl_matrix_complex_get(matrix, i, j)) << " + " << GSL_IMAG(gsl_matrix_complex_get(matrix, i, j)) << "i   "; 
      std::cout << std::scientific << GSL_REAL(gsl_matrix_complex_get(matrix, i, j)) << " + " << GSL_IMAG(gsl_matrix_complex_get(matrix, i, j)) << "i   "; 
    }
    std::cout << std::endl;
  }
}







namespace nusquids {

class nuSQUIDSLIV: public nuSQUIDS {

  private:


    // The parameters of the LIV model
    // int LIV_n;
    // double LIV_cft;

    

    // squids::SU_vector LIVP;
    // std::vector<squids::SU_vector> LIVP_evol;

    squids::SU_vector LIVP_Eindep;
    std::vector<squids::SU_vector> LIVP_Eindep_evol;
    squids::SU_vector LIVP_Edep;
    std::vector<squids::SU_vector> LIVP_Edep_evol;


    void AddToPreDerive(double x){
      for(unsigned int ei = 0; ei < ne; ei++){
        // asumming same mass hamiltonian for neutrinos/antineutrinos
        squids::SU_vector h0 = H0( E_range[ei], 0 );
        // LIVP_evol[ei] = LIVP.Evolve( h0, (x-Get_t_initial()) );

        LIVP_Eindep_evol[ei] = LIVP_Eindep.Evolve( h0, (x-Get_t_initial()) );
        LIVP_Edep_evol[ei] = LIVP_Edep.Evolve( h0, (x-Get_t_initial()) );

      }
    }


    squids:: SU_vector HI(unsigned int ie,unsigned int irho) const {

      squids::SU_vector potential = nuSQUIDS::HI(ie, irho);
      double sign = 1;
      // if (NT==antineutrino){                                      // maybe add (irho == 1 and NT==both)
      //     // antineutrino matter potential flips sign
      //     sign*=(-1);

      // ================= HERE WE ADD THE NEW PHYSICS ===================
      potential += sign * ( E_range[ie] * LIVP_Edep_evol[ie] + LIVP_Eindep_evol[ie]);         //  H_eff_LIV
      // ================= HERE WE ADD THE NEW PHYSICS ===================

      // std::cout << "H_eff : " << potential << std::endl;
      
      return potential;
    }



  void init() {
    /*
      Common init function, called by various constructors
    */

    // Init LIV params to null values
    // LIV_n = 0;
    // LIV_cft = 0;

    // Allocate some matrices
    // LIVP_evol.resize(ne);
    LIVP_Eindep_evol.resize(ne);
    LIVP_Edep_evol.resize(ne);

    for(unsigned int ei = 0; ei < ne; ei++){
      // LIVP_evol[ei] = squids::SU_vector(nsun);
      LIVP_Eindep_evol[ei] = squids::SU_vector(nsun);
      LIVP_Edep_evol[ei] = squids::SU_vector(nsun);
    }

  }


  public:


    //
    // Constructors
    //

    //Default void constructor (never used but required by boost python, specifically RegisterBasicNuSQuIDSPythonBindings, not sure why).
    nuSQUIDSLIV() {}

    /// Multiple energy mode constructor. Basically just passing down to base class constructor and calling init.
    nuSQUIDSLIV(marray<double,1> E_vector,unsigned int numneu, NeutrinoType NT = both, bool iinteraction = false, std::shared_ptr<CrossSectionLibrary> ncs = nullptr)
      : nuSQUIDS(E_vector, numneu, NT, iinteraction, ncs)
      // : nuSQUIDS(E_vector, numneu, NT, false, nullptr)
    { init(); }

    /// Single energy mode constructor. Basically just passing down to base class constructor and calling init.
    nuSQUIDSLIV(unsigned int numneu, NeutrinoType NT = neutrino)
      : nuSQUIDS(numneu, NT)
    { init(); }

    // Load from serialized file (not yet implemented, but required to exist for pybindings)
    nuSQUIDSLIV(std::string hdf5_filename, std::string grp = "/", std::shared_ptr<InteractionStructure> int_struct = nullptr)
      // : nuSQUIDS(hdf5_filename, grp, int_struct)
    {
      throw std::runtime_error("nuSQUIDSLIV::Error::HDF5 serialization not implemented");
    }



    //
    // SME param getters/setters
    //

     void Set_LIVCoefficient(const marray<double,3>& a_mat,const marray<double,3>& c_mat,  const marray<double,3>& e_mat,  double ra_rad, double dec_rad){
      
      // 
      // Assign LIV parameters:
      //
      // dmat is a 3D array of shape (3,3,3) containing a 3x3 matrix for each direction (x,y,z)
      // loop over directions and assign matrices to gsl matrices

      gsl_matrix_complex* a_eV_x = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* a_eV_y = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* a_eV_z = gsl_matrix_complex_calloc(3, 3);

      gsl_matrix_complex* c_tx = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* c_ty = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* c_tz = gsl_matrix_complex_calloc(3, 3);
  

      for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) {
          gsl_matrix_complex_set(a_eV_x, i, j, gsl_complex_rect(a_mat[0][i][j], 0));
          gsl_matrix_complex_set(a_eV_y, i, j, gsl_complex_rect(a_mat[1][i][j], 0));
          gsl_matrix_complex_set(a_eV_z, i, j, gsl_complex_rect(a_mat[2][i][j], 0));
          gsl_matrix_complex_set(c_tx, i, j, gsl_complex_rect(c_mat[0][i][j], 0));
          gsl_matrix_complex_set(c_ty, i, j, gsl_complex_rect(c_mat[1][i][j], 0));
          gsl_matrix_complex_set(c_tz, i, j, gsl_complex_rect(c_mat[2][i][j], 0));
        }}


      // std::cout << "a_eV_y : " << std::endl;
      // print_gsl_matrix(a_eV_y);

      //
      // Define Coordinate system
      //

      // celestial colatitude and longitude
      double theta = M_PI/2 + dec_rad;
      double phi = ra_rad;

      // r vector
      double NX = sin(theta) * cos(phi);
      double NY = sin(theta) * sin(phi);
      double NZ = - cos(theta);



      //
      // Calculate LIV effective Hamiltonian with GSL matrix operations
      //

      // Amplitude equations: (cos(omega_sid L) amplitudes)
      // Ac0 = -NX * ax - NY * ay; 
      // Ac1 = 2 * NX * cxt + 2 * NY * cyt;
      // Const = NZ * az;

      // Declare Ac0, Ac1, and Const matrices
      gsl_matrix_complex* Ac0 = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* Ac1 = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* Const = gsl_matrix_complex_calloc(3, 3);

      // Calculate Ac0, Ac1, and Const matrices
      gsl_matrix_complex_scale(a_eV_x, gsl_complex_rect(-NX, 0));
      gsl_matrix_complex_scale(a_eV_y, gsl_complex_rect(-NY, 0));
      gsl_matrix_complex_memcpy(Ac0, a_eV_x);
      gsl_matrix_complex_add(Ac0, a_eV_y);

      gsl_matrix_complex_scale(c_tx, gsl_complex_rect(2 * NX, 0));
      gsl_matrix_complex_scale(c_ty, gsl_complex_rect(2 * NY, 0));
      gsl_matrix_complex_memcpy(Ac1, c_tx);
      gsl_matrix_complex_add(Ac1, c_ty);

      gsl_matrix_complex_memcpy(Const, a_eV_z);
      gsl_matrix_complex_scale(Const, gsl_complex_rect(NZ, 0));




      // Declare LIVP_Edep and LIVP_Eindep matrices
      gsl_matrix_complex* LIVP_Edep_GSL = gsl_matrix_complex_calloc(3, 3);
      gsl_matrix_complex* LIVP_Eindep_GSL = gsl_matrix_complex_calloc(3, 3);

      // Calculate LIVP_Edep matrix
      gsl_matrix_complex_memcpy(LIVP_Edep_GSL, Ac1);

      // Calculate LIVP_Eindep matrix
      gsl_matrix_complex_memcpy(LIVP_Eindep_GSL, Ac0);
      gsl_matrix_complex_add(LIVP_Eindep_GSL, Const);

      

      // std::cout << "Energy-dependent H_eff : " << std::endl;
      // print_gsl_matrix(LIVP_Edep_GSL);
      // std::cout << "Energy-independent H_eff : " << std::endl;
      // print_gsl_matrix(LIVP_Eindep_GSL);



      // Heff =  Ac0 + Const + E * Ac1 = LIVP_Eindep + E * LIVP_Edep (see HI function above)
      LIVP_Edep = squids::SU_vector(LIVP_Edep_GSL);           // E dependent part of LIVP 
      LIVP_Eindep = squids::SU_vector(LIVP_Eindep_GSL);   // E independent part of LIVP






      
      // // construct SU_vectors from matrices
      // squids::SU_vector ax = squids::SU_vector(a_eV_x);
      // squids::SU_vector ay = squids::SU_vector(a_eV_y);
      // squids::SU_vector az = squids::SU_vector(a_eV_z);
      // squids::SU_vector cxt = squids::SU_vector(c_tx);
      // squids::SU_vector cyt = squids::SU_vector(c_ty);
      // squids::SU_vector czt = squids::SU_vector(c_tz);


      // // celestial colatitude and longitude
      // double theta = M_PI/2 + dec_rad;
      // double phi = ra_rad;

      // // r vector
      // double NX = - sin(theta) * cos(phi);
      // double NY = - sin(theta) * sin(phi);
      // double NZ = - cos(theta);



      // // Amplitude to be multiplied with cos(omega_sid L)
      // squids::SU_vector Ac0 = -NX * ax - NY * ay; 
      // squids::SU_vector Ac1 = 2 * NX * cxt + 2 * NY * cyt;
      // squids::SU_vector Const = NZ * az;


      // // Heff =  Ac0 + Const + E * Ac1 = LIVP_Eindep + E * LIVP_Edep (see HI function above)
      // LIVP_Edep = squids::SU_vector(Ac1);           // E dependent part of LIVP 
      // LIVP_Eindep = squids::SU_vector(Ac0+Const);   // E independent part of LIVP

      // // std::cout << "Energy-dependent H_eff : " << LIVP_Edep << std::endl;
      // // std::cout << "Energy-independent H_eff : " << LIVP_Eindep << std::endl;


       // free allocated matrix
      gsl_matrix_complex_free(a_eV_x);
      gsl_matrix_complex_free(a_eV_y);
      gsl_matrix_complex_free(a_eV_z);
      gsl_matrix_complex_free(c_tx);
      gsl_matrix_complex_free(c_ty);
      gsl_matrix_complex_free(c_tz);
      gsl_matrix_complex_free(Ac0);
      gsl_matrix_complex_free(Ac1);
      gsl_matrix_complex_free(Const);
      gsl_matrix_complex_free(LIVP_Edep_GSL);
      gsl_matrix_complex_free(LIVP_Eindep_GSL);
    }
};



/*
  nuSQUIDSAtm extended to include LIV
*/

class nuSQUIDSLIVAtm : public nuSQUIDSAtm<nuSQUIDSLIV> {

  public:

    // Use the base class constructors
    using nuSQUIDSAtm<nuSQUIDSLIV>::nuSQUIDSAtm;


    // Wrap all the getters/setters

    void Set_LIVCoefficient(const marray<double,3>& a_mat,const marray<double,3>& c_mat,  const marray<double,3>& e_mat,  double ra_rad, double dec_rad){
      for(nuSQUIDSLIV& nsq : this->GetnuSQuIDS()) nsq.Set_LIVCoefficient(a_mat,c_mat,e_mat,ra_rad,dec_rad);
    } 

    void Set_Body(std::shared_ptr<Body> body){
      for(nuSQUIDSLIV& nsq : this->GetnuSQuIDS()) nsq.Set_Body(body);
    }

};



}; // close nusquids namespace

#endif //nusquidslv_h

