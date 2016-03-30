/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize model parameters
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Initialize model parameters
 *****************************************************************************/
void
vic_init(void)
{
    extern all_vars_struct    *all_vars;
    extern size_t              current;
    extern domain_struct       global_domain;
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern filenames_struct    filenames;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;
    extern lake_con_struct    *lake_con;
    extern parameters_struct   param;

    bool                       found;
    bool                       no_overstory;
    char                       locstr[MAXSTRING];
    double                     mean;
    double                     sum;
    double                    *Cv_sum = NULL;
    double                    *dvar = NULL;
    int                       *ivar = NULL;
    size_t                     i;
    size_t                     j;
    size_t                     k;
    size_t                     m;
    size_t                     nveg;
    size_t                     max_numnod;
    int                        veg_class;
    int                        vidx;
    size_t                     d2count[2];
    size_t                     d2start[2];
    size_t                     d3count[3];
    size_t                     d3start[3];
    size_t                     d4count[4];
    size_t                     d4start[4];
    int                        tmp_lake_idx;

    // allocate memory for Cv_sum
    Cv_sum = malloc(local_domain.ncells_active * sizeof(*Cv_sum));
    if (Cv_sum == NULL) {
        log_err("Memory allocation error in vic_init().");
    }

    // allocate memory for variables to be read
    dvar = malloc(local_domain.ncells_active * sizeof(*dvar));
    if (dvar == NULL) {
        log_err("Memory allocation error in vic_init().");
    }
    ivar = malloc(local_domain.ncells_active * sizeof(*ivar));
    if (ivar == NULL) {
        log_err("Memory allocation error in vic_init().");
    }

    // The method used to convert the NetCDF fields to VIC structures for
    // individual grid cells is to read a 2D slice and then loop over the
    // domain cells to assign the values to the VIC structures

    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain.n_ny;
    d2count[1] = global_domain.n_nx;

    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = 1;
    d3count[1] = global_domain.n_ny;
    d3count[2] = global_domain.n_nx;

    d4start[0] = 0;
    d4start[1] = 0;
    d4start[2] = 0;
    d4start[3] = 0;
    d4count[0] = 1;
    d4count[1] = 1;
    d4count[2] = global_domain.n_ny;
    d4count[3] = global_domain.n_nx;

    // start the clock
    current = 0;

    // read_veglib()

    // TBD: Check that options.NVEGTYPES is the right number. Bare soil is
    // the complicating factor here
    for (i = 0; i < local_domain.ncells_active; i++) {
        Cv_sum[i] = 0.;

        for (j = 0; j < options.NVEGTYPES; j++) {
            veg_lib[i][j].NVegLibTypes = options.NVEGTYPES;
            veg_lib[i][j].veg_class = (int) j;
        }
    }

    // overstory
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_int(filenames.veglib, "overstory",
                                 d3start, d3count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_lib[i][j].overstory = ivar[i];
        }
    }

    // rarc
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.veglib, "rarc",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_lib[i][j].rarc = (double) dvar[i];
        }
    }

    // rmin
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.veglib, "rmin",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_lib[i][j].rmin = (double) dvar[i];
        }
    }

    // wind height
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.veglib, "wind_h",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_lib[i][j].wind_h = (double) dvar[i];
        }
    }

    // RGL
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.veglib, "RGL",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_lib[i][j].RGL = (double)dvar[i];
        }
    }

    // rad_atten
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.veglib, "rad_atten",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_lib[i][j].rad_atten = (double) dvar[i];
        }
    }

    // wind_atten
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.veglib, "wind_atten",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_lib[i][j].wind_atten = (double) dvar[i];
        }
    }

    // trunk_ratio
    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.veglib, "trunk_ratio",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_lib[i][j].trunk_ratio = (double) dvar[i];
        }
    }

    // LAI and Wdmax
    if (options.LAI_SRC == FROM_VEGLIB || options.LAI_SRC == FROM_VEGPARAM) {
        for (j = 0; j < options.NVEGTYPES; j++) {
            d4start[0] = j;
            for (k = 0; k < MONTHS_PER_YEAR; k++) {
                d4start[1] = k;
                get_scatter_nc_field_double(filenames.veglib, "LAI",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    veg_lib[i][j].LAI[k] = (double) dvar[i];
                    veg_lib[i][j].Wdmax[k] = param.VEG_LAI_WATER_FACTOR *
                                             veg_lib[i][j].LAI[k];
                }
            }
        }
    }

    // albedo
    if (options.ALB_SRC == FROM_VEGLIB || options.ALB_SRC == FROM_VEGPARAM) {
        for (j = 0; j < options.NVEGTYPES; j++) {
            d4start[0] = j;
            for (k = 0; k < MONTHS_PER_YEAR; k++) {
                d4start[1] = k;
                get_scatter_nc_field_double(filenames.veglib, "albedo",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    veg_lib[i][j].albedo[k] = (double) dvar[i];
                }
            }
        }
    }

    // veg_rough
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < MONTHS_PER_YEAR; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.veglib, "veg_rough",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                veg_lib[i][j].roughness[k] = (double) dvar[i];
            }
        }
    }

    // displacement
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < MONTHS_PER_YEAR; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.veglib, "displacement",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                veg_lib[i][j].displacement[k] = (double) dvar[i];
            }
        }
    }

    // default value for fcanopy
    for (j = 0; j < options.NVEGTYPES; j++) {
        if (options.FCAN_SRC == FROM_DEFAULT) {
            for (k = 0; k < MONTHS_PER_YEAR; k++) {
                for (i = 0; i < local_domain.ncells_active; i++) {
                    veg_lib[i][j].fcanopy[k] = 1.0;
                }
            }
        }
        else if (options.FCAN_SRC == FROM_VEGLIB ||
                 options.FCAN_SRC == FROM_VEGPARAM) {
            d4start[0] = j;
            for (k = 0; k < MONTHS_PER_YEAR; k++) {
                d4start[1] = k;
                get_scatter_nc_field_double(filenames.veglib, "fcanopy",
                                            d4start, d4count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    veg_lib[i][j].fcanopy[k] = (double) dvar[i];
                }
            }
        }
    }

    // read carbon cycle parameters
    if (options.CARBON) {
        // Ctype
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_int(filenames.veglib, "Ctype",
                                     d3start, d3count, ivar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                veg_lib[i][j].Ctype = ivar[i];
                if (veg_lib[i][j].Ctype != PHOTO_C3 &&
                    veg_lib[i][j].Ctype != PHOTO_C4) {
                    log_err("cell %zu veg %zu: Ctype is %d but "
                            "must be either %d (C3) or %d (C4).",
                            i, j, veg_lib[i][j].Ctype,
                            PHOTO_C3, PHOTO_C4);
                }
            }
        }
        // MaxCarboxRate
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.veglib, "MaxCarboxRate",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                veg_lib[i][j].MaxCarboxRate = (double) dvar[i];
                if (veg_lib[i][j].MaxCarboxRate < 0) {
                    log_err("cell %zu veg %zu: MaxCarboxRate is %f "
                            "but must be >= 0.",
                            i, j, veg_lib[i][j].MaxCarboxRate);
                }
            }
        }
        // MaxETransport or CO2Specificity
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.veglib, "MaxiE_or_CO2Spec",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                if (dvar[i] < 0) {
                    log_err("cell %zu veg %zu: MaxE_of_CO2Spec is %f "
                            "but must be >= 0.", i, j, dvar[i]);
                }
                if (veg_lib[i][j].Ctype == PHOTO_C3) {
                    veg_lib[i][j].MaxCarboxRate = (double) dvar[i];
                    veg_lib[i][j].CO2Specificity = 0;
                }
                else if (veg_lib[i][j].Ctype == PHOTO_C4) {
                    veg_lib[i][j].MaxCarboxRate = 0;
                    veg_lib[i][j].CO2Specificity = (double) dvar[i];
                }
            }
        }
        // LightUseEff
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.veglib, "LUE",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                veg_lib[i][j].LightUseEff = (double) dvar[i];
                if (veg_lib[i][j].LightUseEff < 0 ||
                    veg_lib[i][j].LightUseEff > 1) {
                    log_err("cell %zu veg %zu: LightUseEff is %f "
                            "but must be between 0 and 1.",
                            i, j, veg_lib[i][j].LightUseEff);
                }
            }
        }
        // Nscale flag
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_int(filenames.veglib, "Nscale",
                                     d3start, d3count, ivar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                veg_lib[i][j].NscaleFlag = ivar[i];
                if (veg_lib[i][j].NscaleFlag != 0 &&
                    veg_lib[i][j].NscaleFlag != 1) {
                    log_err("cell %zu veg %zu: NscaleFlag is %d but "
                            "must be either 0 or 1.",
                            i, j, veg_lib[i][j].NscaleFlag);
                }
            }
        }
        // Wnpp_inhib
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.veglib, "Wnpp_inhib",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                veg_lib[i][j].Wnpp_inhib = (double) dvar[i];
                if (veg_lib[i][j].Wnpp_inhib < 0 ||
                    veg_lib[i][j].Wnpp_inhib > 1) {
                    log_err("cell %zu veg %zu: Wnpp_inhib is %f "
                            "but must be between 0 and 1.",
                            i, j, veg_lib[i][j].Wnpp_inhib);
                }
            }
        }
        // NPPfactor_sat
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.veglib, "NPPfactor_sat",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                veg_lib[i][j].NPPfactor_sat = (double) dvar[i];
                if (veg_lib[i][j].NPPfactor_sat < 0 ||
                    veg_lib[i][j].NPPfactor_sat > 1) {
                    log_err("cell %zu veg %zu: NPPfactor_sat is %f "
                            "but must be between 0 and 1.",
                            i, j, veg_lib[i][j].NPPfactor_sat);
                }
            }
        }
    }

    // read_soilparam()

    // Nlayer
    options.Nlayer = get_nc_dimension(filenames.soil, "nlayer");

    // Validate Nlayer
    if ((options.FULL_ENERGY ||
         options.FROZEN_SOIL) && options.Nlayer < 3) {
        log_err("You must define at least 3 soil moisture layers to run "
                "the model in FULL_ENERGY or FROZEN_SOIL modes.  "
                "Currently Nlayers is set to  %zu.", options.Nlayer);
    }
    if ((!options.FULL_ENERGY &&
         !options.FROZEN_SOIL) && options.Nlayer < 1) {
        log_err("You must define at least 1 soil moisture layer to run "
                "the model.  Currently Nlayers is set to %zu.",
                options.Nlayer);
    }
    if (options.Nlayer > MAX_LAYERS) {
        log_err("Global file wants more soil moisture layers (%zu) than "
                "are defined by MAX_LAYERS (%d).  Edit vic_driver_shared.h and "
                "recompile.", options.Nlayer, MAX_LAYERS);
    }

    // latitude and longitude
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].lat = local_domain.locations[i].latitude;
        soil_con[i].lng = local_domain.locations[i].longitude;
    }

    // b_infilt
    get_scatter_nc_field_double(filenames.soil, "infilt",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].b_infilt = (double) dvar[i];
    }

    // Ds
    get_scatter_nc_field_double(filenames.soil, "Ds",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].Ds = (double) dvar[i];
    }

    // Dsmax
    get_scatter_nc_field_double(filenames.soil, "Dsmax",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].Dsmax = (double) dvar[i];
    }

    // Ws
    get_scatter_nc_field_double(filenames.soil, "Ws",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].Ws = (double) dvar[i];
    }

    // c
    get_scatter_nc_field_double(filenames.soil, "c",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].c = (double) dvar[i];
    }

    // expt: unsaturated hydraulic conductivity exponent for each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "expt",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].expt[j] = (double) dvar[i];
        }
    }

    // Ksat: saturated hydraulic conductivity for each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "Ksat",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].Ksat[j] = (double) dvar[i];
        }
    }

    // init_moist: initial soil moisture for cold start
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "init_moist",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].init_moist[j] = (double) dvar[i];
        }
    }

    // phi_s
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "phi_s",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].phi_s[j] = (double) dvar[i];
        }
    }

    // elevation: mean grid cell elevation
    get_scatter_nc_field_double(filenames.soil, "elev",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].elevation = (double) dvar[i];
    }

    // depth: thickness for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "depth",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].depth[j] = (double) dvar[i];
        }
    }

    // avg_temp: mean grid temperature
    get_scatter_nc_field_double(filenames.soil, "avg_T",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].avg_temp = (double) dvar[i];
    }

    // dp: damping depth
    get_scatter_nc_field_double(filenames.soil, "dp",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].dp = (double) dvar[i];
    }

    // bubble: bubbling pressure for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "bubble",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].bubble[j] = (double) dvar[i];
        }
    }

    // quartz: quartz content for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "quartz",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].quartz[j] = (double) dvar[i];
        }
    }

    // bulk_dens_min: mineral bulk density for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "bulk_density",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].bulk_dens_min[j] = (double) dvar[i];
        }
    }

    // soil_dens_min: mineral soil density for each soil layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "soil_density",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].soil_dens_min[j] = (double) dvar[i];
        }
    }


    // organic soils
    if (options.ORGANIC_FRACT) {
        // organic
        for (j = 0; j < options.Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.soil, "organic",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                soil_con[i].organic[j] = (double) dvar[i];
            }
        }

        // bulk_dens_org: organic bulk density for each soil layer
        for (j = 0; j < options.Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.soil, "bulk_density_org",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                soil_con[i].bulk_dens_org[j] = (double) dvar[i];
            }
        }

        // soil_dens_org: organic soil density for each soil layer
        for (j = 0; j < options.Nlayer; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.soil, "soil_density_org",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                soil_con[i].soil_dens_org[j] = (double) dvar[i];
            }
        }
    }

    // Wcr: critical point for each layer
    // Note this value is  multiplied with the maximum moisture in each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "Wcr_FRACT",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].Wcr[j] = (double) dvar[i];
        }
    }

    // Wpwp: wilting point for each layer
    // Note this value is  multiplied with the maximum moisture in each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "Wpwp_FRACT",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].Wpwp[j] = (double) dvar[i];
        }
    }

    // rough: soil roughness
    get_scatter_nc_field_double(filenames.soil, "rough",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].rough = (double) dvar[i];
    }

    // snow_rough: snow roughness
    get_scatter_nc_field_double(filenames.soil, "snow_rough",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].snow_rough = (double) dvar[i];
    }

    // annual_prec: annual precipitation
    get_scatter_nc_field_double(filenames.soil, "annual_prec",
                                d2start, d2count, dvar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].annual_prec = (double) dvar[i];
    }

    // resid_moist: residual moisture content for each layer
    for (j = 0; j < options.Nlayer; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.soil, "resid_moist",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].resid_moist[j] = (double) dvar[i];
        }
    }

    // fs_active: frozen soil active flag
    get_scatter_nc_field_int(filenames.soil, "fs_active",
                             d2start, d2count, ivar);
    for (i = 0; i < local_domain.ncells_active; i++) {
        soil_con[i].FS_ACTIVE = (char) ivar[i];
    }

    // spatial snow
    if (options.SPATIAL_SNOW) {
        // max_snow_distrib_slope
        get_scatter_nc_field_double(filenames.soil, "max_snow_distrib_slope",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].max_snow_distrib_slope = (double) dvar[i];
        }

        // frost_slope: slope of frozen soil distribution
        get_scatter_nc_field_double(filenames.soil, "frost_slope",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].frost_slope = (double) dvar[i];
        }
    }

    if (options.COMPUTE_TREELINE) {
        // avgJulyAirTemp: average July air temperature
        get_scatter_nc_field_double(filenames.soil, "avgJulyAirTemp",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].avgJulyAirTemp = (double) dvar[i];
        }
    }

    // Additional processing of the soil variables
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j < options.Nlayer; j++) {
            // compute layer properties
            soil_con[i].bulk_density[j] =
                (1 - soil_con[i].organic[j]) * soil_con[i].bulk_dens_min[j] +
                soil_con[i].organic[j] * soil_con[i].bulk_dens_org[j];
            soil_con[i].soil_density[j] =
                (1 - soil_con[i].organic[j]) * soil_con[i].soil_dens_min[j] +
                soil_con[i].organic[j] * soil_con[i].soil_dens_org[j];
            if (soil_con[i].resid_moist[j] == MISSING) {
                soil_con[i].resid_moist[j] = param.SOIL_RESID_MOIST;
            }
            soil_con[i].porosity[j] = 1 - soil_con[i].bulk_density[j] /
                                      soil_con[i].soil_density[j];
            soil_con[i].max_moist[j] = soil_con[i].depth[j] *
                                       soil_con[i].porosity[j] * MM_PER_M;
            // check layer thicknesses
            if (soil_con[i].depth[j] < MINSOILDEPTH) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_err("Model will not function with layer %zd depth %f < %f "
                        "m.\n%s", j, soil_con[i].depth[j], MINSOILDEPTH,
                        locstr);
            }
        }
        // check relative thickness of top two layers
        if (soil_con[i].depth[0] > soil_con[i].depth[1]) {
            sprint_location(locstr, &(local_domain.locations[i]));
            log_err("Model will not function with layer 0 thicker than layer "
                    "(%f m > %f m).\n%s", soil_con[i].depth[0],
                    soil_con[i].depth[1], locstr);
        }
        // compute maximum infiltration for upper layers
        if (options.Nlayer == 2) {
            soil_con[i].max_infil = (1.0 + soil_con[i].b_infilt) *
                                    soil_con[i].max_moist[0];
        }
        else {
            soil_con[i].max_infil = (1.0 + soil_con[i].b_infilt) *
                                    (soil_con[i].max_moist[0] +
                                     soil_con[i].max_moist[1]);
        }

        // compute soil layer critical and wilting point moisture contents
        for (j = 0; j < options.Nlayer; j++) {
            soil_con[i].Wcr[j] *= soil_con[i].max_moist[j];
            soil_con[i].Wpwp[j] *= soil_con[i].max_moist[j];
            if (soil_con[i].Wpwp[j] > soil_con[i].Wcr[j]) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_err("Calculated wilting point moisture (%f mm) is "
                        "greater than calculated critical point moisture "
                        "(%f mm) for layer %zd."
                        "\n\tIn the soil parameter file, "
                        "Wpwp_FRACT MUST be <= Wcr_FRACT.\n%s",
                        soil_con[i].Wpwp[j], soil_con[i].Wcr[j], j, locstr);
            }
            if (soil_con[i].Wpwp[j] < soil_con[i].resid_moist[j] *
                soil_con[i].depth[j] * MM_PER_M) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_err("Calculated wilting point moisture (%f mm) is "
                        "less than calculated residual moisture (%f mm) for "
                        "layer %zd.\n\tIn the soil parameter file, "
                        "Wpwp_FRACT MUST be >= resid_moist / "
                        "(1.0 - bulk_density/soil_density).\n%s",
                        soil_con[i].Wpwp[j], soil_con[i].resid_moist[j] *
                        soil_con[i].depth[j] * MM_PER_M, j, locstr);
            }
        }

        // validate spatial snow/frost params
        if (options.SPATIAL_SNOW) {
            if (soil_con[i].max_snow_distrib_slope < 0.0) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_err("max_snow_distrib_slope (%f) must be "
                        "positive.\n%s", soil_con[i].max_snow_distrib_slope,
                        locstr);
            }
        }

        if (options.SPATIAL_FROST) {
            if (soil_con[i].frost_slope < 0.0) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_err("frost_slope (%f) must be positive.\n%s",
                        soil_con[i].frost_slope, locstr);
            }
        }

        // If BASEFLOW = NIJSSEN2001 then convert NIJSSEN2001
        // parameters d1, d2, d3, and d4 to ARNO baseflow
        // parameters Ds, Dsmax, Ws, and c
        if (options.BASEFLOW == NIJSSEN2001) {
            j = options.Nlayer - 1;
            soil_con[i].Dsmax = soil_con[i].Dsmax *
                                pow((double)(1. /
                                             (soil_con[i].max_moist[j] -
                                              soil_con[i].Ws)),
                                    -soil_con[i].c) +
                                soil_con[i].Ds * soil_con[i].max_moist[j];
            soil_con[i].Ds = soil_con[i].Ds *
                             soil_con[i].Ws / soil_con[i].Dsmax;
            soil_con[i].Ws = soil_con[i].Ws / soil_con[i].max_moist[j];
        }

        soil_moisture_from_water_table(&(soil_con[i]), options.Nlayer);

        if (options.CARBON) {
            // TBD Remove hardcoded parameter values
            soil_con[i].AlbedoPar = 0.92 * param.ALBEDO_BARE_SOIL - 0.015;
            if (soil_con[i].AlbedoPar < param.PHOTO_ALBSOIPARMIN) {
                soil_con[i].AlbedoPar = param.PHOTO_ALBSOIPARMIN;
            }
        }
    }

    // read_snowband()
    if (options.SNOW_BAND == 1) {
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].AreaFract[0] = 1.;
            soil_con[i].BandElev[0] = soil_con[i].elevation;
            soil_con[i].Pfactor[0] = 1.;
            soil_con[i].Tfactor[0] = 0.;
        }
    }
    else {
        // AreaFract: fraction of grid cell in each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.snowband, "AreaFract",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                soil_con[i].AreaFract[j] = (double) dvar[i];
            }
        }
        // elevation: elevation of each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.snowband, "elevation",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                soil_con[i].BandElev[j] = (double) dvar[i];
            }
        }
        // Pfactor: precipitation multiplier for each snow band
        for (j = 0; j < options.SNOW_BAND; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.snowband, "Pfactor",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                soil_con[i].Pfactor[j] = (double) dvar[i];
            }
        }
        // Run some checks and corrections for soil
        for (i = 0; i < local_domain.ncells_active; i++) {
            // Make sure area fractions are positive and add to 1
            sum = 0.;
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] < 0) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("Negative snow band area fraction (%f) read from "
                            "file\n%s", soil_con[i].AreaFract[j], locstr);
                }
                sum += soil_con[i].AreaFract[j];
            }
            // TBD: Need better check for equal to 1.
            if (sum != 1.) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_warn("Sum of the snow band area fractions does not equal "
                         "1 (%f), dividing each fraction by the sum\n%s",
                         sum, locstr);
                for (j = 0; j < options.SNOW_BAND; j++) {
                    soil_con[i].AreaFract[j] /= sum;
                }
            }
            // check that the mean elevation from the snow bands matches the
            // grid cell mean elevation. If not reset mean
            mean = 0.;
            for (j = 0; j < options.SNOW_BAND; j++) {
                mean += soil_con[i].BandElev[j];
            }
            mean /= options.SNOW_BAND;
            if (fabs(soil_con[i].elevation - soil_con[i].BandElev[j]) > 1.0) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_warn("average band elevation %f not equal to grid_cell "
                         "average elevation %f; setting grid cell elevation "
                         "to average band elevation.\n%s",
                         mean, soil_con[i].elevation, locstr);
                soil_con[i].elevation = (double)mean;
            }
            // Tfactor: calculate the temperature factor
            for (j = 0; j < options.SNOW_BAND; j++) {
                // TBD: Ensure that Tlapse is implemented consistently
                soil_con[i].Tfactor[j] = (soil_con[i].BandElev[j] -
                                          soil_con[i].elevation) *
                                         param.LAPSE_RATE;
            }
            // Pfactor: calculate Pfactor from the precipitation fraction read
            // from file
            // TBD: Ensure that netCDF variable is appropriately named
            sum = 0.;
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].Pfactor[j] < 0.) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("Snow band precipitation fraction (%f) "
                            "must be between 0 and 1.\n%s",
                            soil_con[i].Pfactor[j], locstr);
                }
                if (soil_con[i].Pfactor[j] > 0. &&
                    soil_con[i].AreaFract[j] == 0) {
                    // TBD: Check to make sure whether this check is actually
                    // needed
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("Snow band precipitation fraction (%f) "
                            "should be 0 when the area fraction is "
                            "0. (band = %zd).\n%s",
                            soil_con[i].AreaFract[j], j, locstr);
                }
                sum += soil_con[i].Pfactor[j];
            }
            // TBD: Need better check for equal to 1.
            if (sum != 1.) {
                sprint_location(locstr, &(local_domain.locations[i]));
                log_warn("Sum of the snow band precipitation fractions does "
                         "not equal 1 (%f), dividing each fraction by the "
                         "sum\n%s", sum, locstr);
                for (j = 0; j < options.SNOW_BAND; j++) {
                    soil_con[i].Pfactor[j] /= sum;
                }
            }
            // Pfactor: convert precipitation fraction to Pfactor
            for (j = 0; j < options.SNOW_BAND; j++) {
                if (soil_con[i].AreaFract[j] > 0) {
                    soil_con[i].Pfactor[j] /= soil_con[i].AreaFract[j];
                }
                else {
                    soil_con[i].Pfactor[j] = 0.;
                }
            }
        }
    }

    // logic from compute_treeline()
    for (i = 0; i < local_domain.ncells_active; i++) {
        for (j = 0; j < options.SNOW_BAND; j++) {
            // Lapse average annual July air temperature
            if (soil_con[i].avgJulyAirTemp + soil_con[i].Tfactor[j] <=
                param.TREELINE_TEMPERATURE) {
                // Snow band is above treeline
                soil_con[i].AboveTreeLine[j] = true;
            }
            else {
                soil_con[i].AboveTreeLine[j] = false;
            }
        }
    }

    // read_vegparam()

    // reading the vegetation parameters is slightly more complicated because
    // VIC allocates memory for veg_con only if the vegetation type exists in
    // the grid cell. The veg_con_map_struct is used to provide some of this
    // mapping

    // number of vegetation types - in vic this is defined without the bare soil
    // and the vegetation above the treeline
    for (i = 0; i < local_domain.ncells_active; i++) {
        nveg = veg_con_map[i].nv_active - 1;
        if (options.AboveTreelineVeg >= 0) {
            nveg -= 1;
        }
        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            veg_con[i][j].vegetat_type_num = (int) nveg;
        }
    }

    // Cv: for each vegetation type, read the cover fraction into the mapping
    // structure. Then assign only the ones with a fraction greater than 0 to
    // the veg_con structure

    for (j = 0; j < options.NVEGTYPES; j++) {
        d3start[0] = j;
        get_scatter_nc_field_double(filenames.veg, "Cv",
                                    d3start, d3count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            veg_con_map[i].Cv[j] = (double) dvar[i];
        }
    }

    // do the mapping
    for (i = 0; i < local_domain.ncells_active; i++) {
        k = 0;
        for (j = 0; j < options.NVEGTYPES; j++) {
            if (veg_con_map[i].Cv[j] > 0) {
                veg_con_map[i].vidx[j] = k;
                veg_con[i][k].Cv = veg_con_map[i].Cv[j];
                veg_con[i][k].veg_class = j;
                for (m = 0; m < MONTHS_PER_YEAR; m++) {
                    if (options.LAI_SRC == FROM_VEGLIB ||
                        options.LAI_SRC == FROM_VEGPARAM) {
                        veg_con[i][k].LAI[m] = veg_lib[i][j].LAI[m];
                        veg_con[i][k].Wdmax[m] = param.VEG_LAI_WATER_FACTOR *
                                                 veg_con[i][k].LAI[m];
                    }
                    if (options.ALB_SRC == FROM_VEGLIB ||
                        options.ALB_SRC == FROM_VEGPARAM) {
                        veg_con[i][k].albedo[m] = veg_lib[i][j].albedo[m];
                    }
                    if (options.FCAN_SRC == FROM_VEGLIB ||
                        options.FCAN_SRC == FROM_VEGPARAM) {
                        veg_con[i][k].fcanopy[m] = veg_lib[i][j].fcanopy[m];
                    }
                }
                k++;
            }
            else {
                veg_con_map[i].vidx[j] = NODATA_VEG;
            }
        }
    }

    // zone_depth: root zone depths
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < options.ROOT_ZONES; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.veg, "root_depth",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                vidx = veg_con_map[i].vidx[j];
                if (vidx != NODATA_VEG) {
                    veg_con[i][vidx].zone_depth[k] = (double) dvar[i];
                }
            }
        }
    }

    // zone_fract: root fractions
    for (j = 0; j < options.NVEGTYPES; j++) {
        d4start[0] = j;
        for (k = 0; k < options.ROOT_ZONES; k++) {
            d4start[1] = k;
            get_scatter_nc_field_double(filenames.veg, "root_fract",
                                        d4start, d4count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                vidx = veg_con_map[i].vidx[j];
                if (vidx != NODATA_VEG) {
                    veg_con[i][vidx].zone_fract[k] = (double) dvar[i];
                }
            }
        }
    }

    // calculate root fractions
    for (i = 0; i < local_domain.ncells_active; i++) {
        calc_root_fractions(veg_con[i], &(soil_con[i]));
    }

    // Run some checks and corrections for vegetation
    for (i = 0; i < local_domain.ncells_active; i++) {
        no_overstory = false;
        // Only run to options.NVEGTYPES - 1, since bare soil is the last type
        for (j = 0; j < options.NVEGTYPES - 1; j++) {
            vidx = veg_con_map[i].vidx[j];
            if (vidx != NODATA_VEG) {
                sum = 0;
                for (k = 0; k < options.ROOT_ZONES; k++) {
                    sum += veg_con[i][vidx].zone_depth[k];
                }
                if (sum <= 0) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("Root zone depths must sum to a value greater "
                            "than 0 (sum = %.2f) - Type: %zd.\n%s", sum, j,
                            locstr);
                }
                sum = 0;
                for (k = 0; k < options.ROOT_ZONES; k++) {
                    sum += veg_con[i][vidx].zone_fract[k];
                }
                // TBD: Need better test for not equal to 1.
                if (sum != 1.) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_warn("Root zone fractions sum to more than 1 (%f), "
                             "normalizing fractions.  If the sum is large, "
                             "check your vegetation parameter file.\n%s",
                             sum, locstr);
                    for (k = 0; k < options.ROOT_ZONES; k++) {
                        veg_con[i][vidx].zone_fract[k] /= sum;
                    }
                }
                // check that the vegetation type is defined in the vegetation
                // library
                found = false;
                for (k = 0; k < options.NVEGTYPES; k++) {
                    if (veg_con[i][vidx].veg_class == veg_lib[i][k].veg_class) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("The vegetation class id %i in vegetation tile %i "
                            "from cell %zd is not defined in the vegetation "
                            "library\n%s", veg_con[i][vidx].veg_class, vidx, i,
                            locstr);
                }
                Cv_sum[i] += veg_con[i][vidx].Cv;

                // check for overstory
                if (!veg_lib[i][j].overstory) {
                    no_overstory = true;
                }
            }
        }
        // handle the bare soil portion of the tile
        vidx = veg_con_map[i].vidx[options.NVEGTYPES - 1];
        Cv_sum[i] += veg_con[i][vidx].Cv;

        // handle the vegetation for the treeline option. This is somewhat
        // confusingly handled in VIC. If I am not mistaken, in VIC classic
        // this is handled in the following way:
        //
        // The treeline option is only active if there is more than one snow
        // band and options.COMPUTE_TREELINE is explicitly set in the global
        // file. If the treeline option is active, then there a few cases:
        //
        // 1. The grid cell contains one or more vegetation types that
        // do not have an overstory (either bare soil or vegetation). Nothing
        // further needs to be done to the input. For the elevation bands above
        // the treeline, the values from vegetation with an overstory are simply
        // ignored and the understory and bare ground values are scaled so they
        // cover the entire band. This scaling is done in put_data()
        //
        // 2. The grid cell contains only vegetation with an overstory.
        // In that case a small area of bare soil or vegetation without an
        // overstory must be created.  This will have almost no effect
        // on the results for most elevation bands, but above the treeline, the
        // elevation band will consists entirely of bare soil or the understory
        // vegetation (because of the scaling in put_data(). There are two
        // cases:
        //
        // 2.a. options.AboveTreelineVeg < 0. In that case a small amount of
        // bare soil is created (fraction is 0.001).
        //
        // 2.b. options.AboveTreelineVeg > 0. In that case a small amount of
        // the new vegetation is created (fraction is 0.001). This vegetation
        // should not have an overstory.
        //
        // The tricky parts are:
        //
        // Ensure that the correct number of vegetation types are reflected
        // for each cell.
        //
        // Ensure that bare soil remains the last vegetation type (the one with
        // the highest number). This will seem odd, but that is how it is
        // handled within VIC.
        //
        // Only case 2 needs to be handled explicitly

        if (options.SNOW_BAND > 1 && options.COMPUTE_TREELINE &&
            !no_overstory && Cv_sum[i] == 1.) {
            // Use bare soil above treeline
            if (options.AboveTreelineVeg < 0) {
                for (j = 0; j < options.NVEGTYPES; j++) {
                    vidx = veg_con_map[i].vidx[j];
                    if (vidx != NODATA_VEG) {
                        veg_con[i][vidx].Cv -=
                            0.001 / veg_con[i][vidx].vegetat_type_num;
                    }
                }
                Cv_sum[i] -= 0.001;
            }
            // Use defined vegetation type above treeline
            else {
                for (j = 0; j < options.NVEGTYPES; j++) {
                    vidx = veg_con_map[i].vidx[j];
                    if (vidx != NODATA_VEG) {
                        veg_con[i][vidx].Cv -=
                            0.001 / veg_con[i][vidx].vegetat_type_num;
                        veg_con[i][vidx].vegetat_type_num += 1;
                    }
                }
                veg_con[i][options.NVEGTYPES - 1].Cv = 0.001;
                veg_con[i][options.NVEGTYPES - 1].veg_class =
                    options.AboveTreelineVeg;
                veg_con[i][options.NVEGTYPES - 1].vegetat_type_num =
                    veg_con[i][0].vegetat_type_num;
                // Since root zones are not defined they are copied from another
                // vegetation type.
                for (j = 0; j < options.ROOT_ZONES; j++) {
                    veg_con[i][options.NVEGTYPES - 1].zone_depth[j] =
                        veg_con[i][0].zone_depth[j];
                    veg_con[i][options.NVEGTYPES - 1].zone_fract[j] =
                        veg_con[i][0].zone_fract[j];
                }
                // redo the mapping to ensure that the veg type is active
                k = 0;
                for (j = 0; j < options.NVEGTYPES; j++) {
                    if (veg_con_map[i].Cv[j] > 0) {
                        veg_con_map[i].vidx[j] = k;
                        veg_con[i][k].Cv = veg_con_map[i].Cv[j];
                        veg_con[i][k].veg_class = j;
                        k++;
                    }
                    else {
                        veg_con_map[i].vidx[j] = NODATA_VEG;
                    }
                }
                // check that the vegetation type is defined in the vegetation
                // library
                found = false;
                for (k = 0; k < options.NVEGTYPES; k++) {
                    if (veg_con[i][vidx].veg_class == veg_lib[i][k].veg_class) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("The vegetation class id %i in vegetation tile %i "
                            "from cell %zd is not defined in the vegetation "
                            "library\n%s", veg_con[i][vidx].veg_class, vidx, i,
                            locstr);
                }
                // make sure it has no overstory
                veg_class = veg_con[i][options.NVEGTYPES - 1].veg_class;
                if (veg_lib[i][veg_class].overstory) {
                    sprint_location(locstr, &(local_domain.locations[i]));
                    log_err("Vegetation class %i is defined to have overstory, "
                            "so it cannot be used as the default vegetation "
                            "type for above canopy snow bands.\n%s", veg_class,
                            locstr);
                }
            }
        }

        // If the sum of the tile fractions is not within a tolerance, throw an error
        if (!assert_close_double(Cv_sum[i], 1., 0., 0.001)) {
            sprint_location(locstr, &(local_domain.locations[i]));
            log_err("Cv !=  1.0 (%f) at grid cell %zd. Exiting ...\n%s",
                    Cv_sum[i], i, locstr);
        }
    }

    // read blowing snow parameters
    if (options.BLOWING) {
        // sigma_slope
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.veg, "sigma_slope",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                vidx = veg_con_map[i].vidx[j];
                if (vidx != NODATA_VEG) {
                    veg_con[i][vidx].sigma_slope = (double) dvar[i];
                    if (veg_con[i][vidx].sigma_slope <= 0) {
                        log_err("cell %zu veg %d: deviation of terrain slope "
                                "(sigma_slope) is %f but must be > 0.",
                                i, vidx, veg_con[i][vidx].sigma_slope);
                    }
                }
            }
        }
        // lag_one
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.veg, "lag_one",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                vidx = veg_con_map[i].vidx[j];
                if (vidx != NODATA_VEG) {
                    veg_con[i][vidx].lag_one = (double) dvar[i];
                    if (veg_con[i][vidx].lag_one <= 0) {
                        log_err("cell %zu veg %d: lag_one is %f but "
                                "must be > 0.",
                                i, vidx, veg_con[i][vidx].lag_one);
                    }
                }
            }
        }
        // fetch
        for (j = 0; j < options.NVEGTYPES; j++) {
            d3start[0] = j;
            get_scatter_nc_field_double(filenames.veg, "fetch",
                                        d3start, d3count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                vidx = veg_con_map[i].vidx[j];
                if (vidx != NODATA_VEG) {
                    veg_con[i][vidx].fetch = (double) dvar[i];
                    if (veg_con[i][vidx].fetch <= 1) {
                        log_err("cell %zu veg %d: fetch is %f but "
                                "must be > 1.",
                                i, vidx, veg_con[i][vidx].fetch);
                    }
                }
            }
        }
    }


    // read_lake parameters
    if (options.LAKES) {
        // lake_idx
        get_scatter_nc_field_int(filenames.lakeparam, "lake_idx",
                                 d2start, d2count, ivar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lake_con[i].lake_idx = ivar[i];
            if (!(lake_con[i].lake_idx >= -1 &&
                  lake_con[i].lake_idx <
                  (int) veg_con[i][0].vegetat_type_num)) {
                log_err("cell %zu lake_idx is %d but we must have -1 "
                        "<= lake_idx < Nveg (%zu).\n", i, lake_con[i].lake_idx,
                        veg_con[i][0].vegetat_type_num);
            }
            if (lake_con[i].lake_idx != -1) {
                veg_con[i][lake_con[i].lake_idx].LAKE = 1;
            }
        }

        // numnod
        get_scatter_nc_field_int(filenames.lakeparam, "numnod",
                                 d2start, d2count, ivar);
        max_numnod = 0;
        for (i = 0; i < local_domain.ncells_active; i++) {
            lake_con[i].numnod = (size_t) ivar[i];
            if (lake_con[i].lake_idx == -1) {
                if (lake_con[i].numnod != 0) {
                    log_err("cell %zu lake_idx is %d (lake not present) "
                            "which requires numnod to be 0, but numnod is "
                            "%zu\n", i, lake_con[i].lake_idx,
                            lake_con[i].numnod);
                }
            }
            else if (!(lake_con[i].numnod > 0 &&
                       lake_con[i].numnod < MAX_LAKE_NODES)) {
                log_err("cell %zu numnod is %zu but we must have 1 "
                        "<= numnod < %d.\n", i, lake_con[i].numnod,
                        MAX_LAKE_NODES);
            }
            else if (!(lake_con[i].numnod <= options.NLAKENODES)) {
                log_err("cell %zu numnod is %zu but this exceeds "
                        "the file lake_node dimension length of %zu.\n",
                        i, lake_con[i].numnod, options.NLAKENODES);
            }
            if (lake_con[i].numnod > max_numnod) {
                max_numnod = lake_con[i].numnod;
            }
        }

        // mindepth (minimum depth for which channel outflow occurs)
        get_scatter_nc_field_double(filenames.lakeparam, "mindepth",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lake_con[i].mindepth = (double) dvar[i];
            if (lake_con[i].lake_idx == -1) {
                if (lake_con[i].mindepth != 0) {
                    log_err("cell %zu lake_idx is %d (lake not present) "
                            "which requires mindepth to be 0, but mindepth "
                            "is %f\n", i, lake_con[i].lake_idx,
                            lake_con[i].mindepth);
                }
            }
            else if (lake_con[i].lake_idx != -1 &&
                     !(lake_con[i].mindepth >= 0)) {
                log_err("cell %zu mindepth is %f but must be >= 0.\n",
                        i, lake_con[i].mindepth);
            }
        }

        // wfrac
        get_scatter_nc_field_double(filenames.lakeparam, "wfrac",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lake_con[i].wfrac = (double) dvar[i];
            if (lake_con[i].lake_idx == -1) {
                if (lake_con[i].wfrac != 0) {
                    log_err("cell %zu lake_idx is %d (lake not present) "
                            "which requires wfrac to be 0, but wfrac is "
                            "%f\n", i, lake_con[i].lake_idx, lake_con[i].wfrac);
                }
            }
            else if (lake_con[i].lake_idx != -1 &&
                     !(lake_con[i].wfrac >= 0 && lake_con[i].wfrac <= 1)) {
                log_err("cell %zu wfrac is %f but we must have "
                        "0 <= wfrac <= 1.\n", i, lake_con[i].wfrac);
            }
        }

        // depth_in (initial depth for a cold start)
        get_scatter_nc_field_double(filenames.lakeparam, "depth_in",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lake_con[i].depth_in = (double) dvar[i];
            if (lake_con[i].lake_idx == -1) {
                if (lake_con[i].depth_in != 0) {
                    log_err("cell %zu lake_idx is %d (lake not present) "
                            "which requires depth_in to be 0, but depth_in is "
                            "%f\n", i, lake_con[i].lake_idx,
                            lake_con[i].depth_in);
                }
            }
            else if (lake_con[i].lake_idx != -1 &&
                     !(lake_con[i].depth_in >= 0)) {
                log_err("cell %zu depth_in is %f but must be >= 0.\n",
                        i, lake_con[i].depth_in);
            }
        }

        // rpercent
        get_scatter_nc_field_double(filenames.lakeparam, "rpercent",
                                    d2start, d2count, dvar);
        for (i = 0; i < local_domain.ncells_active; i++) {
            lake_con[i].rpercent = (double) dvar[i];
            if (lake_con[i].lake_idx == -1) {
                if (lake_con[i].rpercent != 0) {
                    log_err("cell %zu lake_idx is %d (lake not present) "
                            "which requires rpercent to be 0, but rpercent is "
                            "%f\n", i, lake_con[i].lake_idx,
                            lake_con[i].rpercent);
                }
            }
            else if (lake_con[i].lake_idx != -1 &&
                     !(lake_con[i].rpercent >= 0 &&
                       lake_con[i].rpercent <= 1)) {
                log_err("cell %zu rpercent is %f but we must have "
                        "0 <= rpercent <= 1.\n", i, lake_con[i].rpercent);
            }
        }

        // lake depth-area relationship
        for (i = 0; i < local_domain.ncells_active; i++) {
            for (j = 0; j <= MAX_LAKE_NODES; j++) {
                lake_con[i].z[j] = 0;
                lake_con[i].Cl[j] = 0;
            }
        }
        if (options.LAKE_PROFILE) {
            for (j = 0; j < max_numnod; j++) {
                d3start[0] = j;

                // basin_depth
                get_scatter_nc_field_double(filenames.lakeparam, "basin_depth",
                                            d3start, d3count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    lake_con[i].z[j] = (double) dvar[i];
                }

                // basin_area
                get_scatter_nc_field_double(filenames.lakeparam, "basin_area",
                                            d3start, d3count, dvar);
                for (i = 0; i < local_domain.ncells_active; i++) {
                    lake_con[i].Cl[j] = (double) dvar[i];
                }
            }
        }
        else {
            // basin_depth
            get_scatter_nc_field_double(filenames.lakeparam, "basin_depth",
                                        d2start, d2count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lake_con[i].z[0] = (double) dvar[i];
            }

            // basin_area
            get_scatter_nc_field_double(filenames.lakeparam, "basin_area",
                                        d2start, d2count, dvar);
            for (i = 0; i < local_domain.ncells_active; i++) {
                lake_con[i].Cl[0] = (double) dvar[i];
            }
        }

        // validate depth-area relationship
        for (i = 0; i < local_domain.ncells_active; i++) {
            // validate top node
            if (lake_con[i].lake_idx == -1) {
                if (lake_con[i].z[0] > 0) {
                    log_err("cell %zu lake_idx is %d (lake not present) "
                            "which requires max depth to be 0, but max depth "
                            "is %f\n", i, lake_con[i].lake_idx,
                            lake_con[i].z[0]);
                }
                if (lake_con[i].Cl[0] > 0) {
                    log_err("cell %zu lake_idx is %d (lake not present) "
                            "which requires max area fraction to be 0, but "
                            "max area fraction is %f\n", i,
                            lake_con[i].lake_idx, lake_con[i].Cl[0]);
                }
            }
            else {
                if (!(lake_con[i].z[0] > 0)) {
                    log_err("cell %zu lake basin max depth is %f but must "
                            "be > 0.\n", i, lake_con[i].z[0]);
                }
                else if (!(lake_con[i].mindepth <= lake_con[i].z[0])) {
                    log_err("cell %zu lake basin mindepth is %f but "
                            "must be <= max depth of %f\n",
                            i, lake_con[i].mindepth, lake_con[i].z[0]);
                }
                if (!(lake_con[i].Cl[0] > 0 && lake_con[i].Cl[0] <= 1)) {
                    log_err("cell %zu lake basin max area fraction is %f but "
                            "we must have 0 < max area fraction < 1\n", i,
                            lake_con[i].Cl[0]);
                }
                if (fabs(1 - lake_con[i].Cl[0] /
                         veg_con[i][lake_con[i].lake_idx].Cv) > 0.01) {
                    log_err("cell %zu lake basin max area fraction is %f but "
                            "must == area fraction of veg tile containing "
                            "lake (%f)\n", i, lake_con[i].Cl[0],
                            veg_con[i][lake_con[i].lake_idx].Cv);
                }
                else {
                    lake_con[i].Cl[0] = veg_con[i][lake_con[i].lake_idx].Cv;
                }
            }

            // valdate other nodes
            if (options.LAKE_PROFILE) {
                for (j = 1; j < lake_con[i].numnod; j++) {
                    if (!(lake_con[i].z[j] > 0 &&
                          lake_con[i].z[j] < lake_con[i].z[j - 1])) {
                        log_err("cell %zu lake basin node %zu depth is %f "
                                "but must be > 0 and < node %zu depth %f\n",
                                i, j, lake_con[i].z[j], j - 1,
                                lake_con[i].z[j - 1]);
                    }
                    if (!(lake_con[i].Cl[j] > 0 &&
                          lake_con[i].Cl[j] < lake_con[i].Cl[j - 1])) {
                        log_err("cell %zu lake basin node %zu area fraction "
                                "is %f but must be > 0 and < node %zu area "
                                "fraction %f\n", i, j, lake_con[i].Cl[j], j - 1,
                                lake_con[i].Cl[j - 1]);
                    }
                }
            }
        }

        // compute other lake parameters here
        for (i = 0; i < local_domain.ncells_active; i++) {
            soil_con[i].cell_area = local_domain.locations[i].area;
            compute_lake_params(&(lake_con[i]), soil_con[i]);
        }
    }

    // initialize structures with default values
    for (i = 0; i < local_domain.ncells_active; i++) {
        nveg = veg_con[i][0].vegetat_type_num;
        initialize_snow(all_vars[i].snow, nveg);
        initialize_soil(all_vars[i].cell, &(soil_con[i]), nveg);
        initialize_veg(all_vars[i].veg_var, nveg);
        if (options.LAKES) {
            tmp_lake_idx = (int)lake_con[i].lake_idx;
            if (tmp_lake_idx < 0) {
                tmp_lake_idx = 0;
            }
            initialize_lake(&(all_vars[i].lake_var), lake_con[i],
                            &(soil_con[i]),
                            &(all_vars[i].cell[tmp_lake_idx][0]), 0);
        }
        initialize_energy(all_vars[i].energy, nveg);
    }

    // cleanup
    free(dvar);
    free(ivar);
}
