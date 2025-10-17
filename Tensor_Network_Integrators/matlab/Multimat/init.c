
#ifdef MPI
#include "mpi.h"
#endif
#include "minip.h" 
#include "mesh.h" 
#include "init.h" 
#include "eos.h" 

static const double mbar             = 1.0e+12; // 1 Mbar = 1.0e+12 Ba
static const double kbar             = 1.0e+09; // 1 Kbar 

int init(char *filename, 
	  int *dims, int *ncells, int *nbdry, Bdry_Type *btype_lower, Bdry_Type *btype_upper, 
          double *xl_prob, double *xr_prob,
          int *nummat, int **matid_ea_mat, int **solid_num_ea_mat, double **gamma_ea_mat,
	  int    **if_fixed_state_ea_mat, double **rho_fixed_state_ea_mat,
	  double **ei_fixed_state_ea_mat, double **pres_fixed_state_ea_mat,  
	  double *tmax, int *ncycle,
          double *dt_viz, int *ncycle_viz, double *courant, double *dt, char **problemname, int mesh_scaling)
{
    char name_implosion[]   = "implosion"; 
    char name_triplepoint[] = "triplepoint"; 
    char name_penetrator[]  = "penetrator"; 
    char name_accelpiston[] = "accelpiston"; 
    char name_shockformation[] = "shockformation"; 
    char name_blowoff[]   = "blowoff";
    char name_advection[] = "advection"; 
    char name_2shocks[] = "2shocks"; 
    char name_sod[] = "sod"; 
    char name_mach60[] = "mach60"; 
    const int nmat_mx = 4;
    const int nreg_mx = 6;
    char inputtype[128], equalsigne[128], probname[128];
    int dim, i, j, k, slen, nmat, nreg, m, r, reg, nb, dir, nshape; 
    int nx, ny, nz; 
    int reg2matids[nreg_mx];
    int ncell[3]; 
    double hot_pres;      // used only for accelerating piston 
    double vel_move_in;   // used only for the problem of 2 shocks.
    double dx, rho0, ei0, p0, p0_check, u0, rad, thick, h; 
    double *c, xl[3], xr[3];
    double rho_ea_reg[nreg_mx], pres_ea_reg[nreg_mx], ei_ea_reg[nreg_mx];
    double v_ea_reg[nreg_mx][3]; 

    double **vel_ea_reg;
  
    Region_Shape reg_shape[nreg_mx];

    double param_outer[4], param_inner[4]; 
    double param_left[6]; 
    double param_block[6]; 
    double param_hot_reg[16];
    double param_piston[16];
    double param_sphere[4]; 
    double param_hex[24]; // (x0,y0,z0), (x1,y1,z1), ..., (x7,y7,z7)  
    double param_cyl[8];  // (rad1,xc1,yc1,zc1, rad2,xc2,yc2,zc2) 

    double pi, angle, cosa, sina, tana, x0_low, x0_up, x1_low, x1_up; 
    FILE *fp; 

    //! Add mesh_scaling for benchmarking
    int scale = 1 << mesh_scaling;
    // nx *= scale;
    // ny *= scale;
    // nz *= scale;

    fp = fopen(filename, "r");
    if (!filename) { 
	printf("ERROR: file %s doesn't exist.\n", filename);
	return 1; 
    } 
    fscanf(fp, "%s%s%s", inputtype, equalsigne, probname);
    while((inputtype[0] == ' ') || (inputtype[0] == '#') || (inputtype[0] == '\\') || 
          (inputtype[0] == '$') || (inputtype[0] == '@') || (inputtype[0] == '&')) { 
          fscanf(fp, "%s%s%s", inputtype, equalsigne, probname); 
    }  
    printf("%s : %s\n",  inputtype, probname); 
    if (!strcmp(probname, name_accelpiston)) {
        fscanf(fp, "%s%s%lf", inputtype, equalsigne,  &hot_pres);
	printf("%s : %e\n",  inputtype, hot_pres); 
    } 
    if (!strcmp(probname, name_2shocks)) {
        fscanf(fp, "%s%s%lf", inputtype, equalsigne, &vel_move_in);
	printf("%s : %e\n",  inputtype, vel_move_in);
    } 
    fscanf(fp, "%s%s%d", inputtype, equalsigne, &dim);
    printf("%s : %d\n",  inputtype, dim);
    fscanf(fp, "%s%s%d", inputtype, equalsigne, &nx);
    nx = nx*scale;
    printf("%s : %d\n",  inputtype, nx);
    fscanf(fp, "%s%s%d", inputtype, equalsigne, &ny);
    ny = ny * scale;
    printf("%s : %d\n",  inputtype, ny);
    fscanf(fp, "%s%s%d", inputtype, equalsigne, &nz);
    nz = nz * scale;
    printf("%s : %d\n",  inputtype, nz);
 
    fscanf(fp, "%s%s%d", inputtype, equalsigne, &nmat); 
    printf("%s : %d\n",  inputtype, nmat); 
    assert((nmat > 0) && (nmat <= nmat_mx)); 
    fscanf(fp, "%s%s%lf", inputtype, equalsigne,  tmax);  
    printf("%s : %e\n",  inputtype, *tmax); 
    fscanf(fp, "%s%s%d", inputtype, equalsigne, ncycle);
    printf("%s : %d\n",  inputtype, *ncycle); 
    fscanf(fp, "%s%s%lf", inputtype, equalsigne, dt_viz);
    printf("%s : %e\n",  inputtype, *dt_viz);
    fscanf(fp, "%s%s%d", inputtype, equalsigne, ncycle_viz);
    printf("%s : %d\n",  inputtype, *ncycle_viz); 
    fscanf(fp, "%s%s%lf", inputtype, equalsigne, courant);
    printf("%s : %e\n",  inputtype, *courant); 
    fscanf(fp, "%s%s%lf", inputtype, equalsigne, dt);
    printf("%s : %e\n",  inputtype, *dt);



    slen = strlen(probname) + 1;
    *problemname = (char *) malloc(slen * sizeof(char));
    strcpy(*problemname, probname);
 
    nb = NBDRY; 
   
    *gamma_ea_mat    = (double *) malloc(nmat * sizeof(double));
    *matid_ea_mat    = (int *) malloc(nmat * sizeof(int)); 
    *solid_num_ea_mat = (int *) malloc(nmat * sizeof(int));
    *if_fixed_state_ea_mat = (int *) malloc(nmat * sizeof(int)); 

    *rho_fixed_state_ea_mat = (double *) malloc(nmat * sizeof(double));
    *ei_fixed_state_ea_mat  = (double *) malloc(nmat * sizeof(double)); 
    *pres_fixed_state_ea_mat = (double *) malloc(nmat * sizeof(double)); 

    for (m = 0; m < nmat; m++) {
        (*matid_ea_mat)[m]    = m;  
	(*solid_num_ea_mat)[m] = 0;
	(*if_fixed_state_ea_mat)[m]   = 0; 
	(*rho_fixed_state_ea_mat)[m]  = 0.0;
	(*ei_fixed_state_ea_mat)[m]   = 0.0;
	(*pres_fixed_state_ea_mat)[m] = 0.0; 
    }	 
    for (m = 0; m < nmat; m++) {
        (*gamma_ea_mat)[m] = 1.4;
    }
    reg_shape[0].type       = shape_universe;
    reg_shape[0].parameters = NULL;
    
    for (reg = 0; reg < nreg_mx; reg++) {
        for (i = 0; i < 3; i++) {
            v_ea_reg[reg][i] = 0.0;
        }
    }
    if (!strcmp(probname, name_implosion)) { 

//      mat0: hot gas 
//      mat1: metal   
//      mat2: air  

        (*if_fixed_state_ea_mat)[0] = 1;  // hot gas  
//        (*if_fixed_state_ea_mat)[1] = 1;  // metal  

        for (m = 0; m < nmat; m++) {
            (*gamma_ea_mat)[m] = 1.4;  
        }
        nreg = 3;   

        reg2matids[0] = 0;  // hot gas  
        reg2matids[1] = 1;  // metal  
        reg2matids[2] = 2;  // air  
         
//      hot gas 
        rho_ea_reg[0]  = 1.0;     // g/cm^3 
        pres_ea_reg[0] = 5.0 * kbar;
        ei0 = pres_ea_reg[0]/((*gamma_ea_mat)[0] - 1.0);

        (*rho_fixed_state_ea_mat)[0]  = rho_ea_reg[0];
        (*pres_fixed_state_ea_mat)[0] = pres_ea_reg[0];
        (*ei_fixed_state_ea_mat)[0]   = ei0;
        
//      metal 
//      rho0 = 8.96;   // copper  
        rho0 = 4.510;    // titanium  
        p0   = 1.0;
        rho_polynominal(p0, &rho0);
        rho_ea_reg[1]  = rho0;
        pres_ea_reg[1] = p0;

        (*solid_num_ea_mat)[1] = 2;
        e_solid((*solid_num_ea_mat)[1], rho0, p0, &ei0);
//        ei0 = pres_ea_reg[1]/((*gamma_ea_mat)[1] - 1.0);
        (*rho_fixed_state_ea_mat)[1]  = rho0;
        (*pres_fixed_state_ea_mat)[1] = p0;
        (*ei_fixed_state_ea_mat)[1]   = ei0;
        
//      air 
        rho_ea_reg[2]  = 1.3e-03;     // g/cm^3 
        pres_ea_reg[2] = p0;

        xl[0] = 0.0; 
        xr[0] = 10.0;
        dx    = (xr[0] - xl[0])/(double)nx;
        xl[1] = 0.0; 
        xr[1] = xl[1] + dx * (double)ny;
        xl[2] = 0.0;
        xr[2] = xl[2] + dx * (double)nz;
 
        ncell[0] = nx; 
        ncell[1] = ny;
        ncell[2] = nz;

//      outer boundary of metal 
        reg_shape[1].type       = shape_sphere;
        reg_shape[1].parameters = param_outer;
        if (dim == 2) {
            param_outer[0] = 4.0;   // radius 
            param_outer[1] = 5.0;   // ctr_x 
            param_outer[2] = 5.0;   // ctr_y ;
        }
        else if (dim == 3) {
            param_outer[0] = 4.0;  // radius 
            param_outer[1] = 5.0;
            param_outer[2] = 5.0; 
            param_outer[3] = 5.0;
        }
//      inner boundary of metal

        reg_shape[2].type       = shape_sphere;
        reg_shape[2].parameters = param_inner;
        if (dim == 2) {
            param_inner[0] = 3.5;   // radius 
            param_inner[1] = 5.0;   // ctr_x 
            param_inner[2] = 5.0;   // ctr_y 
        }   
        else if (dim == 3) {
            param_inner[0] = 3.0;  // radius 
            param_inner[1] = 5.0;
            param_inner[2] = 5.0;  
            param_inner[3] = 5.0;
        }   
    } 
    else if (!strcmp(probname, name_triplepoint)) { 

        double param_upper[4], param_lower[4];

	xl[0] = 0.0;
        xr[0] = 7.0;
        ncell[0] = 1400;
        xl[1] = 0.0;
        xr[1] = 3.0;
        ncell[1] = 600;

	param_lower[0] = 1.0;
        param_lower[1] = 0.0;
        param_lower[2] = 7.0;
        param_lower[3] = 1.5;

	param_upper[0] = 1.0;
	param_upper[1] = 1.5;
	param_upper[2] = 7.0;
	param_upper[3] = 3.0;

	reg_shape[1].type = shape_rectangular;
	reg_shape[1].parameters = param_lower;

    reg_shape[2].type = shape_rectangular;
    reg_shape[2].parameters = param_upper;

	nmat = 2;
	(*gamma_ea_mat)[0] = 1.5;
	(*gamma_ea_mat)[1] = 1.4; 
	
    nreg = 3;
    reg2matids[0] = 0;   // left region 
    reg2matids[1] = 1;   // lower region  
	reg2matids[2] = 0;   // upper region  

//      left region  
	rho_ea_reg[0]  = 1.0;
        pres_ea_reg[0] = 1.0;

//      lower region   
        rho_ea_reg[1]  = 1.0;
        pres_ea_reg[1] = 0.1;

//      upper region 
        rho_ea_reg[2]  = 0.125;
        pres_ea_reg[2] = 0.1;

	for (r = 0; r < nreg; r++) {
             m = reg2matids[r]; 
            ei_ea_reg[r] = pres_ea_reg[r]/((*gamma_ea_mat)[m] - 1.0);
        } 
    } 
    else if (!strcmp(probname, name_penetrator)) {

//      mat0: gas 
//      mat1: copper plate 
//      mat2: copper cylinder 

       nreg = 3;
//     reg0 :  gas
//     reg1 :  plate 
//     reg2 :  cylinder 
    
        xl[0] = -10.0;
        xr[0] =  10.0;
        ncell[0] = 20;   
        xl[1] = 0.0;
        xr[1] = 15.0;
        ncell[1] = 15;  
        xl[2] = 0.0;
        xr[2] = 15.0;
        ncell[2] = 15;  
//                         . 
//                       . . 
//                     .   .   
//                   .     .  
//                 .       .  
//               .         .   
//             .   .       .  
//           . a      .    . 
//         .------------.--.----
//        x0_low           x0_up  
       
	 pi = 2.0 * acos(0.0);
         angle = 60.0 * (pi / 180.0);
         sina = sin(angle);
         cosa = sqrt(1.0 - sina * sina);
         tana = sina / cosa; 
         x0_low = 0.0;
         x0_up  = x0_low+ (xr[2] - xl[2])/tana; 
         thick = 2.0;  // thickness of the plate
         x1_low = x0_low + thick/sina;
         x1_up  = x1_low + (xr[2] - xl[2])/tana;   

         c = param_hex;
         c[0] = x0_low;
         c[1] = xl[1];
         c[2] = xl[2];
         c += dim;
         c[0] = x1_low;
         c[1] = xl[1];
         c[2] = xl[2];
         c += dim;
         c[0] = x1_low;
         c[1] = xr[1];
         c[2] = xl[2]; 
         c += dim;
         c[0] = x0_low;
         c[1] = xr[1];
         c[2] = xl[2];
// 
         c += dim;
         c[0] = x0_up;
         c[1] = xl[1];
         c[2] = xr[2];
         c += dim;
         c[0] = x1_up;
         c[1] = xl[1];
         c[2] = xr[2];
         c += dim;
         c[0] = x1_up;
         c[1] = xr[1];
         c[2] = xr[2];
         c += dim;
         c[0] = x0_up;
         c[1] = xr[1];
         c[2] = xr[2];

//       for cylinder  

	 rad = 4.0; 
	 h   = 10.0;
         c = param_cyl;  
         c[0] = rad;
         c++;
         c[0] = - h;
         c[1] = 0.5 *(xl[1] + xr[1]);
         c[2] = 0.5 *(xl[2] + xr[2]);
         c += 3;
	 c[0] = rad; 
	 c++; 
         c[0] = 0.0;
         c[1] = 0.5 *(xl[1] + xr[1]);
         c[2] = 0.5 *(xl[2] + xr[2]);

	 nreg = 3;
	 reg2matids[0] = 0;  // gas background
         reg2matids[1] = 1;  // copper plate 
         reg2matids[2] = 2;  // copper cylinder  
			     
//       gas background 
         rho_ea_reg[0]  = 1.0e-03;    // g/cc 
         pres_ea_reg[0] = 1.01326e+06; // 1 atm 

//       copper plate 
         (*solid_num_ea_mat)[1] = 2;

         rho0 = 8.96;   // copper  
         p0   = pres_ea_reg[0];
         e_solid((*solid_num_ea_mat)[1], rho0, p0, &ei0);
         p_solid((*solid_num_ea_mat)[1], rho0, ei0, &p0_check);
         assert(fabs(p0_check - p0)/p0 < 1.0e-06);
	 rho_ea_reg[1]  = rho0;
         pres_ea_reg[1] = p0;

//       copper cyliner 

	 rho_ea_reg[2]  = rho0;
         pres_ea_reg[2] = p0;

	 (*solid_num_ea_mat)[2] = 2;

	 reg_shape[1].type       = shape_hex ;
         reg_shape[1].parameters = param_hex; 

	 reg_shape[2].type       = shape_cylinder;
         reg_shape[2].parameters = param_cyl;
  
	 u0 = 1.0e+06; 
	 v_ea_reg[2][0] = u0;
    } 
    else if (!strcmp(probname, name_accelpiston)) { 

//      mat0: vaccum
//      mat1: gas
//      mat2: piston

        for (m = 0; m < nmat; m++) {
            (*gamma_ea_mat)[m] = 1.4;  
        }
        (*solid_num_ea_mat)[2] = 2;

//        (*if_fixed_state_ea_mat)[1] = 1;  // gas  
//        (*if_fixed_state_ea_mat)[2] = 1;  // piston  

//        nreg = 4;   // inslude the sphere  
        nreg = 3;   // exclude sphere   

        reg2matids[0] = 0;  // vocuum background
        reg2matids[1] = 1;  // gas 
        reg2matids[2] = 2;  // pistone  
        reg2matids[3] = 0;  // circle   
         
//      vacuum 
        rho_ea_reg[0]  = 0.0e-00;
        pres_ea_reg[0] = 0.0e+00;

//      gas 
        rho_ea_reg[1]  = 1.0;      // g/cm^3 
        pres_ea_reg[1] = hot_pres; // dyne/cm^2  1 kbar = e+09 dyne/cm^2

	(*rho_fixed_state_ea_mat)[1] = rho_ea_reg[1];
	(*pres_fixed_state_ea_mat)[1] = pres_ea_reg[1];
	(*ei_fixed_state_ea_mat)[1]  = pres_ea_reg[1]/((*gamma_ea_mat)[1] - 1.0);

//      piston 
//        rho0 = 8.96;   // copper  
        rho0 = 4.510;    // titanium  
        p0   = 0.0;

        e_solid((*solid_num_ea_mat)[2], rho0, p0, &ei0);
        p_solid((*solid_num_ea_mat)[2], rho0, ei0, &p0_check); 
       
        rho_ea_reg[2]  = rho0;  
        pres_ea_reg[2] = p0;
	(*rho_fixed_state_ea_mat)[2]  = rho0;
	(*pres_fixed_state_ea_mat)[2] = p0; 
	(*ei_fixed_state_ea_mat)[2]   = ei0;

//      half sphere, vacuum  
        rho_ea_reg[3]  = 0.0;
        pres_ea_reg[3] = 0.0;

        xl[0] = -50.0; 
        xr[0] =  50.0;
        dx    = (xr[0] - xl[0])/(double)nx;
        xl[1] = 0.0; 
        xr[1] = xl[1] + dx * (double)ny;
        xl[2] = 0.0;
        xr[2] = xl[2] + dx * (double)nz;
 
        ncell[0] = nx; 
        ncell[1] = ny;
        ncell[2] = nz;

        reg_shape[1].type       = shape_rectangular;
        reg_shape[1].parameters = param_left;
        if (dim == 2) { 
            param_left[0] = xl[0];      // lower left corner  
            param_left[1] = xl[1];

            param_left[2] = 0.5 *(xl[0] + xr[0]);
            param_left[3] = xr[1];
        }
        else if (dim == 3) {
            param_left[0] = xl[0];
            param_left[1] = xl[1];
            param_left[2] = xl[2];

            param_left[3] = 0.5 *(xl[0] + xr[0]);
            param_left[4] = xr[1];
            param_left[5] = xr[2];
        }
        reg_shape[2].type       = shape_rectangular;
        reg_shape[2].parameters = param_block;
        if (dim == 2) {
            param_block[0] = 0.5 *(xl[0] + xr[0]);      // lower left corner  
            param_block[1] = xl[1];

            param_block[2] = 0.5 *(xl[0] + xr[0]) + 5.0;     // upper right corner 
            param_block[3] = xr[1];
        }
        else if (dim == 3) {
            param_block[0] = 0.5 *(xl[0] + xr[0]);
            param_block[1] = xl[1];
            param_block[2] = xl[2];

            param_block[3] = 0.5 *(xl[0] + xr[0]) + 5.0; 
            param_block[4] = xr[1];
            param_block[5] = xr[2];
        }
       reg_shape[3].type       = shape_sphere;
       reg_shape[3].parameters = param_sphere;

        if (dim == 2) { 
	    param_sphere[0] = 2.5;   // radius 
	    param_sphere[1] = 0.5 *(xl[0] + xr[0]) + 5.0;   // ctr x
            param_sphere[2] = 0.5 *(xl[1] + xr[1]); 
        }
        else if (dim == 3) {
            param_sphere[0] = 2.5;  // radius 
            param_sphere[1] = 0.5 *(xl[0] + xr[0]) + 5.0; 
            param_sphere[2] = 0.5 *(xl[1] + xr[1]);
            param_sphere[3] = 0.5 *(xl[2] + xr[2]);
        }
    } 
    else if (!strcmp(probname, name_shockformation)) { 
  
        for (m = 0; m < nmat; m++) {
            (*gamma_ea_mat)[m] = 1.4;  
        }
	nreg = 2; 
        for (reg = 0; reg < nreg; reg++) { 
            reg2matids[reg] = MIN(reg, nmat-1); 
        } 
        rho_ea_reg[0]  = 1.0e-06;   // g/cm^3 
        pres_ea_reg[0] = 1.0e+04;   // dyne/cm^2,   1 Mbar = e+12 dyne/cm^2  

        rho_ea_reg[1]  = 1.0e+02;   // g/cm"3. heavy eneough 
        pres_ea_reg[1] = 1.0e+04;  

        xl[0] =  0.0; 
        xr[0] = 50.0;
        ncell[0] = 500;
        xl[1] = 0.0; 
        xr[1] = 1.0;
        ncell[1] = 10;  
        xl[2] = 0.0;
        xr[2] = 1.0; 
        ncell[2] = 10; 

        reg_shape[1].type       = shape_rectangular;
        reg_shape[1].parameters = param_left;
        if (dim == 2) { 
            param_left[0] = 0.0;      // lower left corner  
            param_left[1] = 0.0;

            param_left[2] = 0.2;     // upper right corner 
            param_left[3] = 1.0;
        }
        else if (dim == 3) {
            param_left[0] = 0.0;
            param_left[1] = 0.0;
            param_left[2] = 0.0;

            param_left[3] = 0.2;
            param_left[4] = 1.0;
            param_left[5] = 1.0;
        }
        v_ea_reg[1][0] = 1.18e+05;   // 0.118 cm/microsec 
    } 
    else if (!strcmp(probname, name_blowoff)) { 
  
        for (m = 0; m < nmat; m++) {
            (*gamma_ea_mat)[m] = 1.667;  
        }
	nreg = 2; 
        for (reg = 0; reg < nreg; reg++) { 
            reg2matids[reg] = MIN(reg, nmat-1); 
        } 
        rho_ea_reg[0]  = 0.0;
        pres_ea_reg[0] = 0.0;
        
        rho_ea_reg[1]  = 0.20;   // g/cm^3 
        pres_ea_reg[1] = 2.4e+10; // dyne/cm^2  0.024 Mbar  

        xl[0] = -20.0;
        xr[0] =  20.0;
        dx = (xr[0] - xl[0])/(double)nx;
        xl[1] = 0.0;
        xr[1] = xl[1] + dx * (double)ny;
        xl[2] = 0.0;
        xr[2] = xl[2] + dx * (double)nz;

        ncell[0] = nx;
        ncell[1] = ny;
        ncell[2] = nz;

        reg_shape[1].type       = shape_rectangular;
        reg_shape[1].parameters = param_left;
        if (dim == 2) { 
            param_left[0] = xl[0];      // lower left corner  
            param_left[1] = xl[1];

            param_left[2] = 0.5*(xl[0] + xr[0]);     // upper right corner 
            param_left[3] = xr[1];
        }
        else if (dim == 3) {
            param_left[0] = xl[0];
            param_left[1] = xl[1]; 
            param_left[2] = xl[2];

            param_left[3] = 0.5*(xl[0] + xr[0]);     // upper right corner 
            param_left[4] = xr[1];
            param_left[5] = xr[2];
        }
    } 
    else if (!strcmp(probname, name_advection)) { 

        int if_move_to_right = 1;    // 1 for moving to right, 0 for mving to the left 
        int dir_advection    = 0;    // 0, or 1, or 2 for the axis along which the advection happens. 
        if (dir_advection >= dim) dir_advection = 0; 

        (*gamma_ea_mat)[0] = 1.4;
        (*gamma_ea_mat)[1] = 1.6;

        nreg = 2; 
        reg2matids[0] = 0; 
        if (nmat == 1) {  
            reg2matids[1] = 0;
        }
        else { 
            for (reg = 1; reg < nreg; reg++) { 
                m = reg;
                reg2matids[reg] = (*matid_ea_mat)[m];
            }
        } 
        rho_ea_reg[0]  = 1.0;
        pres_ea_reg[0] = 1.0;

        rho_ea_reg[1]  = 2.0;
        pres_ea_reg[1] = 1.0;
 
	    if (dir_advection == 0) {
                xl[0] = 0.0; 
		xr[0] = 1.0;
                dx = (xr[0] - xl[0])/(double)nx; 
		xl[1] = 0.0;
		xl[2] = 0.0;
		xr[1] = xl[1] + ((double) ny) * dx; 
                xr[2] = xl[2] + ((double) nz) * dx;
	    }
	    else if (dir_advection == 1) {
		xl[1] = 0.0;
                xr[1] = 1.0;
                dx = (xr[1] - xl[1])/(double)ny;
		xl[0] = 0.0;
                xl[2] = 0.0;
                xr[0] = xl[0] + ((double) nx) * dx;
                xr[2] = xl[2] + ((double) nz) * dx;
	    }
	    else if (dir_advection == 2) { 
                xl[2] = 0.0;
                xr[2] = 1.0;
                dx = (xr[2] - xl[2])/(double)nz;
		xl[1] = 0.0;
                xl[0] = 0.0;
                xr[1] = xl[1] + ((double) ny) * dx;
                xr[0] = xl[0] + ((double) nx) * dx;
            }
        if (dim == 2) {
	    ncell[0] = nx;
	    ncell[1] = ny;
	    ncell[2] = 1;
	}
	else if (dim == 3) { 
            ncell[0] = nx;
            ncell[1] = ny;
            ncell[2] = nz;
        }
//        reg_shape[1].type       = shape_sphere;
//        reg_shape[1].parameters = param_inner;
        reg_shape[1].type       = shape_rectangular;
        reg_shape[1].parameters = param_left;

        if (dim == 2) { 
           if (dir_advection == 0) {  
                    param_inner[0] = 0.05;   // radius  
                    param_inner[1] = 0.25;   // ctr_x    
		    param_inner[2] = 0.1;    // ctr_y 
					     
                    param_left[0] = 0.0;
                    param_left[1] = 0.0;

                    param_left[2] = 0.5;
                    param_left[3] = 0.2;
            }
            else if (dir_advection == 1) {
		    param_inner[0] = 0.05;   // radius
                    param_inner[1] = 0.1;    // ctr_x 
                    param_inner[2] = 0.25;   // ctr_y 
            }
        }
        else if (dim == 3) { 
            if (dir_advection == 0) {
                param_inner[0] = 0.05;   // radius  
                param_inner[1] = 0.25;   // ctr_x    
                param_inner[2] = 0.1;    // ctr_y 
                param_inner[3] = 0.1;   // ctr_z 	
            }
            else if (dir_advection == 1) {
                param_inner[0] = 0.05;   // radius
                param_inner[1] = 0.1;    // ctr_x 
                param_inner[2] = 0.25;   // ctr_y 
                param_inner[3] = 0.1;    // ctr_z    
            }
            else if (dir_advection == 2) {
		param_inner[0] = 0.05;   // radius
                param_inner[1] = 0.1;    // ctr_x
                param_inner[2] = 0.1;   // ctr_y
                param_inner[3] = 0.25;    // ctr_z
            }
        }  
        if (if_move_to_right) {   
            for (reg = 0; reg < nreg; reg++) { 
                v_ea_reg[reg][dir_advection] = 1.0;
            }
        }
        else { 
            for (reg = 0; reg < nreg; reg++) {
                v_ea_reg[reg][dir_advection] = -1.0;
            } 
        }
    } 
    if (!strcmp(probname, name_2shocks)) { 
        for (m = 0; m < nmat; m++) {
	    (*gamma_ea_mat)[m] = 1.4; 
        }
	(*gamma_ea_mat)[0] = 1.4;
        (*gamma_ea_mat)[1] = 1.6;

        nreg = 2; 
        reg2matids[0] = 0; 
        if (nmat == 1) {  
            reg2matids[1] = 0;
        }
        else { 
            for (reg = 1; reg < nreg; reg++) { 
                m = reg;
                reg2matids[reg] = (*matid_ea_mat)[m];
            }
        } 
        rho_ea_reg[0]  = 1.0;
        pres_ea_reg[0] = 1.0;

        rho_ea_reg[1]  = 1.0;
        pres_ea_reg[1] = 1.0;

        if (dim == 2) {
            xl[0] = -0.5;
	    xr[0] = 0.5;
	    dx = (xr[0] - xl[0])/(double)nx;  
	    xl[1] = 0.0;
            xr[1] = xl[1] + dx * (double)ny;
            ncell[0] = nx;
            ncell[1] = ny;
            ncell[2] = 1;
        } 
        else if (dim == 3) {
            xl[0] = -0.5;
            xr[0] = 0.5;
            dx = (xr[0] - xl[0])/(double)nx;
            xl[1] = 0.0;
            xr[1] = xl[1] + dx * (double)ny;
            xl[2] = 0.0;
            xr[2] = xl[2] + dx * (double)nz;
            ncell[0] = nx;
            ncell[1] = ny;
            ncell[2] = nz;
        }
        reg_shape[1].type       = shape_rectangular;
        reg_shape[1].parameters = param_left;

        if (dim == 2) { 
            param_left[0] = -0.5;    // lower left corner  
            param_left[1] = 0.0;

            param_left[2] = 0.0;     // upper right corner 
            param_left[3] = 0.2;
        }
        else if (dim == 3) { 
            param_left[0] = -0.5;
            param_left[1] = 0.0;
            param_left[2] = 0.0;

            param_left[3] = 0.0;
            param_left[4] = 0.2;
            param_left[5] = 0.2;
        }  
        v_ea_reg[0][0] = - vel_move_in;
        v_ea_reg[1][0] =   vel_move_in; 
    } 
    else if (!strcmp(probname, name_sod)) { 
        for (m = 0; m < nmat; m++) {
	    (*gamma_ea_mat)[m] = 1.4;
        }
        nreg          = 2;
        reg2matids[0] = 0;
        reg2matids[1] = nmat - 1;

        rho_ea_reg[0]  = 0.125;
        rho_ea_reg[1]  = 1.0;
        pres_ea_reg[0] = 0.1;
        pres_ea_reg[1] = 1.0;

        reg_shape[1].type       = shape_rectangular;
        reg_shape[1].parameters = param_left;

        xl[0] = 0.0;
        xr[0] = 1.0;
        dx = (xr[0] - xl[0])/(double)nx;
        xl[1] = 0.0;
        xr[1] = xl[1] + dx * (double)ny;
        xl[2] = 0.0;
        xr[2] = xl[2] + dx * (double)nz;  
        ncell[0] = nx;
        ncell[1] = ny;
        ncell[2] = nz;
    
        if (dim == 2) { 
            param_left[0] = xl[0];
            param_left[1] = xl[0];
	    
	    param_left[2] = 0.5 *(xr[0] + xl[0]);
	    param_left[3] = xr[1];
       }
       else if (dim == 3) {
            param_left[0] = xl[0];
            param_left[1] = xl[1];
            param_left[2] = xl[2];
            
            param_left[3] = 0.5 * (xr[0] + xl[0]);
            param_left[4] = xr[1];
            param_left[5] = xr[2]; 
       }
    } 
    else if (!strcmp(probname, name_mach60)) { 
        (*gamma_ea_mat)[0] = 1.4;
//      (*gamma_ea_mat)[1] = 1.4;

        nreg          = 2;
        reg2matids[0] = 0;
        reg2matids[1] = 0;

        rho_ea_reg[0]  = 1.0;
        pres_ea_reg[0] = 0.1;
	rho_ea_reg[1]  = 3.99666;
        pres_ea_reg[1] = 449.975;

        reg_shape[1].type = shape_rectangular;
        reg_shape[1].parameters = param_left;

        if (dim == 2) {
            xl[0] = 0.0;
            xl[1] = 0.0;
            xr[0] = 1.0;
            xr[1] = 0.2;
            ncell[0] = 100;
            ncell[1] = 20;

            // reg_shape[1].type = shape_rectangular; 
            param_left[0] = 0.0;
            param_left[1] = 0.0;

            param_left[2] = 0.405;
            param_left[3] = 0.2;
       }
       else if (dim == 3) {
            xl[0] = 0.0;
            xl[1] = 0.0;
            xl[2] = 0.0;
            xr[0] = 1.0;
            xr[1] = 0.2;
            xr[2] = 0.2;
            ncell[0] = 100;
            ncell[1] = 20;
            ncell[2] = 20;

            // reg_shape[1].type = shape_rectangular;
            param_left[0] = 0.0;
            param_left[1] = 0.0;
            param_left[2] = 0.0;

            param_left[3] = 0.405;
            param_left[4] = 0.2;
            param_left[5] = 0.2;
       }
       v_ea_reg[0][0] = -18.3661; 
    }
    for (reg = 0; reg < nreg; reg++) {
        m = reg2matids[reg];
	if (!(*solid_num_ea_mat)[m]) { 
            ei_ea_reg[reg] = pres_ea_reg[reg] /((*gamma_ea_mat)[m] - 1.0);
        }
	else { 
            e_solid((*solid_num_ea_mat)[m], rho_ea_reg[reg], pres_ea_reg[reg], ei_ea_reg + reg);  
        } 
    }
    set_mesh(dim, xl, xr, ncell, nb, nmat);

    vel_ea_reg = (double **) malloc(nreg * sizeof(double *));
    for (reg = 0; reg < nreg; reg++) { 
        vel_ea_reg[reg] = v_ea_reg[reg];
    }
    for (dir = 0; dir < dim; dir++) {
        btype_lower[dir] = bdry_transmitted;
        btype_upper[dir] = bdry_transmitted;
    }
    set_mesh_mat(dim, btype_lower, btype_upper, 
		 nmat, *solid_num_ea_mat, *gamma_ea_mat,  
		 nreg, reg2matids, reg_shape, 
                 rho_ea_reg, pres_ea_reg, ei_ea_reg, vel_ea_reg); 
    

                 
    free(vel_ea_reg);

    for (dir = 0; dir < dim; dir++) {
        xl_prob[dir] = xl[dir]; 
        xr_prob[dir] = xr[dir]; 
    } 
    *dims = dim; 
    *nbdry = nb;
    for (dir = 0; dir < dim; dir++) { 
        ncells[dir] = ncell[dir];
    } 
    *nummat = nmat;

    fclose(fp);

    return 0;
} 
