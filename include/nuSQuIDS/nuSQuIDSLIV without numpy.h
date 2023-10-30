#ifndef nusquidslv_H
#define nusquidslv_H

#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <iomanip>
#include <math.h>
#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>



// Printing matrix function
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

    // Sidereal vs isotropic
    bool directional = false;

    // The parameters of the LIV model
    int LIV_n;
    double LIV_cft;


    double ra_rad;
    double dec_rad;
    double sme_a_x;
    double sme_a_y;
    double sme_a_z;
    double sme_c_tx;
    double sme_c_ty;
    double sme_c_tz;

    
    


    //
    //Coordinates
    //

    // celestial colatitude and longitude
    double theta = M_pi/2 + dec_rad;
    double phi = ra 

    
    // spherical coordinates unit propagation vectors
    
    // r vector
    double NX = sin(theta) * cos(phi)
    double NY = sin(theta) * sin(phi)
    double = cos(theta)
    
    // theta vector
    double ThetaX = cos(theta) * cos(phi)
    double ThetaY = cos(theta) * sin(phi)
    double ThetaZ = -sin(theta)  
    
    // phi vector
    double PhiX = -sin(phi)
    double PhiY = cos(phi)
    double PhiZ = 0

    //
    //  Mass independent operators
    //
    
    // Amplitude to be multiplied with sin(omega_sid L)
    double As0 = -NY * a_eV_x + NX * a_eV_y
    double As1 = + 2 * NY * c_tx - 2 * NX * c_ty
    double As = As0 + E * As1

    // Amplitude to be multiplied with cos(omega_sid L)
    double Ac0 = - NX * a_eV_x - NY * a_eV_y
    double Ac1 = 2 * NX * c_tx + 2 * NY * c_ty
    double Ac = Ac0 + E * Ac1



    // double LIV_a_eV;
    // double LIV_c;

//TODO: Fix if-statement for sidereal
    // The LIV potential
    // void def_LIVP() {
    //   if (directional) {
    //     squids::SU_vector LIVPX;
    //     squids::SU_vector LIVPY;
    //     squids::SU_vector LIVPZ;
    //     std::vector<squids::SU_vector> LIVPX_evol;
    //     std::vector<squids::SU_vector> LIVPY_evol;
    //     std::vector<squids::SU_vector> LIVPZ_evol;
    //   } else {
    //     squids::SU_vector LIVP;
    //     std::vector<squids::SU_vector> LIVP_evol;
    //   }
    // }  

    squids::SU_vector LIVP;
    std::vector<squids::SU_vector> LIVP_evol;



    void AddToPreDerive(double x){
      for(unsigned int ei = 0; ei < ne; ei++){
        // asumming same mass hamiltonian for neutrinos/antineutrinos
        squids::SU_vector h0 = H0( E_range[ei], 0 );
        // LIVP_evol[ei] = LIVP.Evolve( h0, (x-Get_t_initial()) );

        // //TODO: Fix if-statement for sidereal
        if (directional) {
          std::cout << "Sidereal not yet implemented. Changed to isotropic" << std::endl;
        //   LIVPX_evol[ei] = LIVPX.Evolve( h0, (x-Get_t_initial()) );
        //   LIVPY_evol[ei] = LIVPY.Evolve( h0, (x-Get_t_initial()) );
        //   LIVPZ_evol[ei] = LIVPZ.Evolve( h0, (x-Get_t_initial()) );
          LIVP_evol[ei] = LIVP.Evolve( h0, (x-Get_t_initial()) );
        } else {
          LIVP_evol[ei] = LIVP.Evolve( h0, (x-Get_t_initial()) );
        }

        
      }
    }


    // squids::SU_vector H0(double E, unsigned int irho) const {
    //   squids::SU_vector h0 = nuSQUIDS::H0( E, irho );
    //   squids::SU_vector h0_liv = h0 + ( LIV_cft * pow(E, LIV_n) );
    //   return h0_liv;
    // }


    squids:: SU_vector HI(unsigned int ie,unsigned int irho) const {
      //std::cout << "+++ HI : start" << std::endl;
      squids::SU_vector potential = nuSQUIDS::HI(ie, irho);
      double sign = 1;
      // if ((irho == 1 and NT==both) or NT==antineutrino){
      //     // antineutrino matter potential flips sign
      //     sign*=(-1);
      // }
      // ================= HERE WE ADD THE NEW PHYSICS ===================
      potential += sign*pow(E_range[ie], LIV_n) * LIVP_evol[ie]; // <- super important line here is where all the physics is set
      // ================= HERE WE ADD THE NEW PHYSICS ===================
      //std::cout << "+++ HI : end" << std::endl;
      return potential;
    }

  void init() {
    /*
      Common init function, called by various constructors
    */

    //std::cout << "+++ nuSQUIDSLIV::init" << std::endl;

    // Init LIV params to null values
    LIV_n = 0;
    LIV_cft = 0;

    // Allocate some matrices
    LIVP_evol.resize(ne);
    for(unsigned int ei = 0; ei < ne; ei++){
      LIVP_evol[ei] = squids::SU_vector(nsun);
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

    void  Set_Directional_SME( //TODO: all the SME parameters should be matrices
      double ra_rad,
      double dec_rad,
      double a_eV_x,
      double a_eV_y,
      double a_eV_z,
      double c_tx,
      double c_ty,
      double c_tz,
      double E){
        
      //
      //Coordinates
      //

      // celestial colatitude and longitude
      double theta = M_pi/2 + dec_rad;
      double phi = ra;

      
      // spherical coordinates unit propagation vectors
      
      // r vector
      double NX = sin(theta) * cos(phi);
      double NY = sin(theta) * sin(phi);
      double NZ = cos(theta);
      

      //
      //  Mass independent operators
      //
      
      // Amplitude to be multiplied with sin(omega_sid L)
      double As0 = -NY * a_eV_x + NX * a_eV_y;
      double As1 = + 2 * NY * c_tx - 2 * NX * c_ty;
      double As = As0 + E * As1;

      // Amplitude to be multiplied with cos(omega_sid L)
      double Ac0 = - NX * a_eV_x - NY * a_eV_y;
      double Ac1 = 2 * NX * c_tx + 2 * NY * c_ty;
      double Ac = Ac0 + E * Ac1;

      double Const = NZ * a_eV_z  
      H_eff = Ac + const 


    }




    void Set_LIVCoefficient(double a_eV, double c_t, std::string field_orientation ){
      
      gsl_complex zero{{ 0.0 , 0.0 }};
      gsl_complex one{{ 1.0 , 0.0 }};
      gsl_complex two{{ 2.0 , 0.0 }};
      gsl_complex three{{ 3.0 , 0.0 }};
      
      gsl_complex a{{ a_eV , 0.0 }}; //Only using real part right now
      gsl_complex c{{ c_t , 0.0 }}; //Only using real part right now


      gsl_matrix_complex * Flavor_structure = gsl_matrix_complex_calloc(3,3);
      gsl_matrix_complex_set(Flavor_structure, 2, 2, one);



      // 
      // Construct LV_Field_structure
      //

      // untill properly implemented:
      field_orientation = "x";

      gsl_matrix_complex * LV_Field_structure = gsl_matrix_complex_calloc(3,3);  //orientation of the LV field

      if (field_orientation == "x"){
        gsl_matrix_complex_set(LV_Field_structure, 0, 0, one);                     //along the x-axis
      }
      else if (field_orientation == "y"){
        gsl_matrix_complex_set(LV_Field_structure, 1, 1, one);                     //along the y-axis
      }
      else if (field_orientation == "z"){
        gsl_matrix_complex_set(LV_Field_structure, 2, 2, one);                     //along the z-axis
      }
      else if (field_orientation == "xy"){
        gsl_matrix_complex_set(LV_Field_structure, 0, 0, one);                     //along the x-axis
        gsl_matrix_complex_set(LV_Field_structure, 1, 1, one);                     //along the y-axis
        gsl_matrix_complex_scale(LV_Field_structure, gsl_complex_inverse(gsl_complex_sqrt(two))); //normalize
      }
      else if (field_orientation == "xz"){
        gsl_matrix_complex_set(LV_Field_structure, 0, 0, one);                     //along the x-axis
        gsl_matrix_complex_set(LV_Field_structure, 2, 2, one);                     //along the z-axis
        gsl_matrix_complex_scale(LV_Field_structure, gsl_complex_inverse(gsl_complex_sqrt(two))); //normalize
      }
      else if (field_orientation == "yz"){
        gsl_matrix_complex_set(LV_Field_structure, 1, 1, one);                     //along the y-axis
        gsl_matrix_complex_set(LV_Field_structure, 2, 2, one);                     //along the z-axis
        gsl_matrix_complex_scale(LV_Field_structure, gsl_complex_inverse(gsl_complex_sqrt(two))); //normalize
      }
      else if (field_orientation == "xyz"){
        gsl_matrix_complex_set(LV_Field_structure, 0, 0, one);                     //along the x-axis
        gsl_matrix_complex_set(LV_Field_structure, 1, 1, one);                     //along the y-axis
        gsl_matrix_complex_set(LV_Field_structure, 2, 2, one);                     //along the z-axis
        gsl_matrix_complex_scale(LV_Field_structure, gsl_complex_inverse(gsl_complex_sqrt(three))); //normalize
      }
      else
        std::cout << "Error: invalid field_orientation" << std::endl;


      gsl_matrix_complex_set(LV_Field_structure, 0, 0, one);                     //along the x-axis

      
      // Define LIV_coef_a matrix and store:  a_eV * (Flavor_structure * LV_Field_structure), no transposition
      gsl_matrix_complex * LIV_coef_a = gsl_matrix_complex_alloc(3, 3);
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, a, Flavor_structure, LV_Field_structure, one, LIV_coef_a);

      // Define LIV_coef_c matrix and store:  c_t * (Flavor_structure * LV_Field_structure), no transposition
      gsl_matrix_complex * LIV_coef_c = gsl_matrix_complex_alloc(3, 3);
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, c, Flavor_structure, LV_Field_structure, one, LIV_coef_c);

      
      


      // gsl_matrix_complex * SME_Coefficients = gsl_matrix_complex_calloc(3,3); //coefficient of the LV field

       // defining a complex matrix M=flavor_matrix*sme_matrix = Heff
       gsl_matrix_complex * M = gsl_matrix_complex_calloc(3,3); //TODO check num nu 
       gsl_matrix_complex_set(M, 2, 2, c);


       LIVP = squids::SU_vector(M);

      // Print heff_matrix
      //  std::cout << "Matrix:" << std::endl;
      //  print_gsl_matrix(M);

       // rotate from flavor to mass basis
       // LIVP.RotateToB1(params);

       // free allocated matrix
       gsl_matrix_complex_free(M);


       // c_params = lv_params;
       // LIV_params_set = true;
    }

    void Set_LIVEnergyDependence(int n){
      LIV_n = n;
    }

};

} // close nusquids namespace

#endif //nusquidslv_h

