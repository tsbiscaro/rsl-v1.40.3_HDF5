/*
HDF5 handling added by:

Thiago S. Biscaro - thiago.biscaro@cptec.inpe.br

Divisao de Satelites e Sistemas Ambientais
Centro de Previsao de Tempo e Estudos Climaticos
Instituto Nacional de Pesquisas Espaciais
Rodovia Presidente Dutra km 40
Cachoeira Paulista/SP 12630000 Brazil

*/


/*
    NASA/TRMM, Code 910.1.q
    This is the TRMM Office Radar Software Library.
    Copyright (C) 1996, 1997
            Mike Kolander
            Space Applications Corporation
            Vienna, Virginia

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public
    License along with this library; if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, US.
*/
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LIBHDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
/* RSL function and structure definitions. */
#include "rsl.h"

#ifdef HAVE_LIBHDF5

#define DEFINED_VOLUMES 9
#define VERSION1_8_0
#define MAX_HDF_STR 255
#define MAXSWEEPS 40

int return_vartype(char *var)
   {
     char *nome_var[DEFINED_VOLUMES] = {"Z", "UZ", "V", "W", "ZDR", "PhiDP", 
                                        "uPhiDP", "KDP", "RhoHV"};
     int vars[DEFINED_VOLUMES] = {CZ_INDEX, DZ_INDEX, VR_INDEX, 
                                  SW_INDEX, DR_INDEX, PH_INDEX,
                                  NP_INDEX, KD_INDEX, RH_INDEX};
   int i;
   
   for (i = 0; i < DEFINED_VOLUMES; i++)
     {
       if (0 == strcmp(var, nome_var[i]))
         {
	   return (vars[i]);
         }
     }
   return (-1);
   }


Radar *RSL_hdf5_ODIM_to_radar(char *infile)
  {
  hid_t file, scan, vol, attr, dataspace, memspace, how,
    header, elev0, elev1, azim0, azim1, what, memtype, where,
    time;
  int status;
  int variavel;
  
  int n_sweep, n_ray, n_bin, n_volume;
  char grp_name[MAX_HDF_STR];
  int rank;
  
  char moments[MAX_RADAR_VOLUMES][MAX_HDF_STR];
  char date_tmp[MAX_HDF_STR];
  char version_str[MAX_HDF_STR];
  
  int pulse;
  double r_step = 0, r_sample = 0, r_start = 0,
     horiz_beam = 0, verti_beam = 0, nyq_vel = 0,
     height = 0, lat = 0, lon = 0, min = 0, sec = 0,
     wavelength = 0, elev = 0, pulsewidth = 0,
     prf = 0, prf2 = 0, gain = 0, undetect = 0,
     space = 0, offset = 0, nodata = 0, time_dbl = 0;
  
  double *azAngleEnd = NULL, *azAngleIni = NULL;

  int ndims = 0;

  float vmax = -32;
  float vmin = 100;
  float vtemp = 0;
  
  float r_min = 0, r_max = 0, scale_factor = 0, elev_header = 0;
  char arquivo[MAX_HDF_STR];
  hsize_t dims[1], sz, width, height2, planes, npals;
  int *data_out = NULL;
  double *ray_elev0 = NULL, *ray_elev1 = NULL;
  double *ray_azim0 = NULL, *ray_azim1 = NULL;
  unsigned long long *timestamp = NULL;
  char date_str[9], time_str[7],
     var_name[8], *data_size = NULL;
  
  unsigned char *image = NULL;
  float *image32 = NULL;
  int i;
  
  int NSWEEPS = 0, NVOLUMES = 0, NRAYS = 0, NBINS = 0;
  
  time_t scan_time;
  H5T_class_t class_id;
  hsize_t dims2;
  size_t sz2;
  

  Radar *radar = NULL;
  Sweep *sweep = NULL;
  Volume *volume[MAX_RADAR_VOLUMES];
  Ray *ray = NULL;
  struct tm *cal_time = NULL;
  
  static float (*f)(Range x);
  static Range (*invf)(float x);

  int voltype = -1, valid_vol = 0;

  int use_volume[MAX_RADAR_VOLUMES];
  
  strcpy(arquivo, infile);

   
  file = H5Fopen(arquivo, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0)
     {
     printf("ERRO\n");
     return (-1);
     }


  /*saves the number of scans (elevations) in this file*/
  for (n_sweep=0; n_sweep < MAXSWEEPS; n_sweep++)
      {
      memset(grp_name, 0, sizeof(grp_name));
      sprintf(grp_name, "dataset%d", n_sweep + 1);
#ifdef VERSION1_8_0
      scan = H5Gopen(file, grp_name, H5P_DEFAULT);
#else
      scan = H5Gopen(file, grp_name);
#endif
      if (scan < 0)
         {
         NSWEEPS = n_sweep;
         break;
         }
      }
     
  
#ifdef VERSION1_8_0   
  how = H5Gopen(file, "/how", H5P_DEFAULT);
  what = H5Gopen(file, "/what", H5P_DEFAULT);
  where = H5Gopen(file, "/where", H5P_DEFAULT);
#else
   how = H5Gopen(file, "/how");
   what = H5Gopen(file, "/what");
   where = H5Gopen(file, "/where");
#endif


   attr = H5Aopen(what, "version", H5P_DEFAULT);
   memtype = H5Aget_type(attr);
   memset(version_str, 0, sizeof(version_str));
   status = H5Aread(attr, memtype, &version_str);
   H5Tclose(memtype);
   H5Aclose(attr);
   
   if (0 != strncmp(version_str, "H5Rad 2.1", 9))
      {
      printf("Versao invalida do HDF (deveria ser 2.1)\n");
      return 0;
      }
   

   attr = H5Aopen(where, "height", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &height);
   H5Aclose(attr);

   attr = H5Aopen(where, "lat", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &lat);
   H5Aclose(attr);

   attr = H5Aopen(where, "lon", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &lon);
   H5Aclose(attr);

   H5Gclose(where);

   /*reading date/time information*/
   attr = H5Aopen(what, "date", H5P_DEFAULT);
   memtype = H5Aget_type(attr);
   memset(date_str, 0, sizeof(date_str));
   status = H5Aread(attr, memtype, &date_str);
   H5Tclose(memtype);
   H5Aclose(attr);

   attr = H5Aopen(what, "time", H5P_DEFAULT);
   memtype = H5Aget_type(attr);
   memset(time_str, 0, sizeof(time_str));
   status = H5Aread(attr, memtype, &time_str);
   H5Tclose(memtype);
   H5Aclose(attr);

   H5Gclose(what);
 
   attr = H5Aopen(how, "beamwH", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &horiz_beam);
   H5Aclose(attr);
   
   attr = H5Aopen(how, "beamwV", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &verti_beam);
   H5Aclose(attr);
   
   H5Gclose(how);
   
   radar = RSL_new_radar(MAX_RADAR_VOLUMES);
   for (n_volume = 0; n_volume < MAX_RADAR_VOLUMES; n_volume++)
      {   
      volume[n_volume] = RSL_new_volume(NSWEEPS);
      }
   
   /*writes lat/lon and height info to the radar header*/
   radar->h.height = height;

   radar->h.latd = (int) lat;
   radar->h.lond = (int) lon;
   lat = fabs(lat);
   lon = fabs(lon);
   radar->h.latm = (int) ((lat - abs(radar->h.latd)) * 60);
   radar->h.lonm = (int) ((lon - abs(radar->h.lond)) * 60);
   radar->h.lats = (int) ((((lat - abs(radar->h.latd)) * 60) -
                           radar->h.latm) * 60);
   radar->h.lons = (int) ((((lon - abs(radar->h.lond)) * 60) -
                           radar->h.lonm) * 60);

   date_tmp[0] = date_str[0];
   date_tmp[1] = date_str[1];
   date_tmp[2] = date_str[2];
   date_tmp[3] = date_str[3];
   date_tmp[4] = '\0';
   radar->h.year = atoi(date_tmp);
   date_tmp[0] = date_str[4];
   date_tmp[1] = date_str[5];
   date_tmp[2] = '\0';
   radar->h.month = atoi(date_tmp);
   date_tmp[0] = date_str[6];
   date_tmp[1] = date_str[7];
   date_tmp[2] = '\0';
   radar->h.day = atoi(date_tmp);
   date_tmp[0] = time_str[0];
   date_tmp[1] = time_str[1];
   date_tmp[2] = '\0';
   radar->h.hour = atoi(date_tmp);
   date_tmp[0] = time_str[2];
   date_tmp[1] = time_str[3];
   date_tmp[2] = '\0';
   radar->h.minute = atoi(date_tmp);
   date_tmp[0] = time_str[4];
   date_tmp[1] = time_str[5];
   date_tmp[2] = '\0';
   radar->h.sec = atoi(date_tmp);

   memset(use_volume, -1, MAX_RADAR_VOLUMES*sizeof(int));
   
   for (n_sweep=0; n_sweep < NSWEEPS; n_sweep++)
      {
      memset(grp_name, 0, sizeof(grp_name));
      sprintf(grp_name, "dataset%d", n_sweep + 1);
#ifdef VERSION1_8_0
      scan = H5Gopen(file, grp_name, H5P_DEFAULT);
#else
      scan = H5Gopen(file, grp_name);
#endif
      if (scan < 0)
         {
         break;
         }
      
#ifdef VERSION1_8_0   
      how = H5Gopen(scan, "how", H5P_DEFAULT);
      what = H5Gopen(scan, "what", H5P_DEFAULT);
      where = H5Gopen(scan, "where", H5P_DEFAULT);
#else
      how = H5Gopen(scan, "how");
      what = H5Gopen(scan, "what");
      where = H5Gopen(scan, "where");
#endif
      
      attr = H5Aopen(where, "elangle", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &elev);
      H5Aclose(attr);
      
      attr = H5Aopen(where, "nbins", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_INT, &NBINS);
      H5Aclose(attr);
      
      attr = H5Aopen(where, "nrays", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_INT, &NRAYS);
      H5Aclose(attr);
      
      attr = H5Aopen(where, "rscale", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &r_step);
      H5Aclose(attr);
      
      attr = H5Aopen(where, "rstart", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &r_start);
      H5Aclose(attr);
      
      H5Gclose(where);
      
      attr = H5Aopen(how, "wavelength", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &wavelength);
      H5Aclose(attr);
      
      attr = H5Aopen(how, "highprf", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &prf);
      H5Aclose(attr);
      
      attr = H5Aopen(how, "lowprf", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &prf2);
      H5Aclose(attr);
      
      attr = H5Aopen(how, "pulsewidth", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &pulsewidth);
      H5Aclose(attr);
      
      attr = H5Aopen(how, "NI", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &nyq_vel);
      H5Aclose(attr);
      
      attr = H5Aopen(how, "beamwH", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &horiz_beam);
      H5Aclose(attr);
      
      attr = H5Aopen(how, "beamwV", H5P_DEFAULT);
      status = H5Aread(attr, H5T_NATIVE_DOUBLE, &verti_beam);
      H5Aclose(attr);
      
      H5Gclose(how);
      
      attr = H5Aopen(what, "startdate", H5P_DEFAULT);
      memtype = H5Aget_type(attr);
      memset(date_str, 0, sizeof(date_str));
      status = H5Aread(attr, memtype, &date_str);
      H5Tclose(memtype);
      H5Aclose(attr);
      
      attr = H5Aopen(what, "starttime", H5P_DEFAULT);
      memtype = H5Aget_type(attr);
      memset(time_str, 0, sizeof(time_str));
      status = H5Aread(attr, memtype, &time_str);
      H5Tclose(memtype);
      H5Aclose(attr);
      
      H5Gclose(what);

      
      NVOLUMES = 0;
      
      for (n_volume = 0; n_volume < MAX_RADAR_VOLUMES; n_volume++)
         {
         memset(grp_name, 0, sizeof(grp_name));
         sprintf(grp_name, "data%d", n_volume + 1);
#ifdef VERSION1_8_0
         vol = H5Gopen(scan, grp_name, H5P_DEFAULT);
#else
         vol = H5Gopen(scan, grp_name);
#endif
	 if (vol < 0)
            {
            break;
            }
         
         volume[n_volume]->sweep[n_sweep] = RSL_new_sweep(NRAYS);
         volume[n_volume]->sweep[n_sweep]->h.sweep_num = n_sweep + 1;
         volume[n_volume]->sweep[n_sweep]->h.beam_width =
            (verti_beam + horiz_beam)/2;
         volume[n_volume]->sweep[n_sweep]->h.vert_half_bw =
            verti_beam / 2;
         volume[n_volume]->sweep[n_sweep]->h.horz_half_bw =
            horiz_beam / 2;
         volume[n_volume]->sweep[n_sweep]->h.nrays = NRAYS;
	 volume[n_volume]->sweep[n_sweep]->h.elev = elev;

#ifdef VERSION1_8_0   
         how = H5Gopen(vol, "how", H5P_DEFAULT);
         what = H5Gopen(vol, "what", H5P_DEFAULT);
#else
         how = H5Gopen(vol, "how");
         what = H5Gopen(vol, "what");
#endif   
            
         attr = H5Aopen(what, "quantity", H5P_DEFAULT);
         memtype = H5Aget_type(attr);
         /*string size in BYTES*/
         memset(var_name, 0, sizeof(var_name));
         status = H5Aread(attr, memtype, &var_name);
         H5Tclose(memtype);
         H5Aclose(attr);


         voltype = return_type_var(var_name);
         
         if (0 > voltype)
            {
            H5Gclose(how);
            H5Gclose(what);
            continue;
            }

         /*marca que esse volume tem uma variavel valida*/
         use_volume[n_volume] = voltype;
         
         attr = H5Aopen(what, "gain", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_DOUBLE, &gain);
         H5Aclose(attr);
            
         attr = H5Aopen(what, "nodata", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_DOUBLE, &nodata);
         H5Aclose(attr);
            
         attr = H5Aopen(what, "offset", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_DOUBLE, &offset);
         H5Aclose(attr);
            
         attr = H5Aopen(what, "undetect", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_DOUBLE, &undetect);
         H5Aclose(attr);
            
         H5Gclose(what);
            
         attr = H5Aopen(how, "startazA", H5P_DEFAULT);
         space = H5Aget_space (attr);
         NRAYS = H5Sget_simple_extent_npoints(space);
         azAngleIni = (double *) calloc(NRAYS, sizeof(double));
         memtype = H5Tcopy(H5T_IEEE_F64LE);
         status = H5Aread (attr, memtype, azAngleIni);
         H5Tclose(memtype);
         H5Sclose (space);
         H5Aclose(attr);
            
         attr = H5Aopen(how, "stopazA", H5P_DEFAULT);
         space = H5Aget_space (attr);
         NRAYS = H5Sget_simple_extent_npoints(space);
         azAngleEnd = (double *) calloc(NRAYS, sizeof(double));
         memtype = H5Tcopy(H5T_IEEE_F64LE);
         status = H5Aread (attr, memtype, azAngleEnd);
         H5Tclose(memtype);
         H5Sclose (space);
         H5Aclose(attr);
         H5Gclose(how);
            
         status = H5LTget_dataset_info(vol, "data", &dims2, &class_id, &sz2);
         if (H5T_INTEGER == class_id)
            {
            image = (unsigned char *) malloc(NRAYS*NBINS*sizeof(unsigned char));
            status = H5LTread_dataset (vol, "data", H5T_NATIVE_UCHAR, image);
            }
         else
            {
            image32 = (float *) malloc(NRAYS*NBINS*sizeof(float));
            status = H5LTread_dataset (vol, "data", H5T_NATIVE_FLOAT, image32);
            }
         
         switch (voltype)
            {
            case DZ_INDEX:
               f = DZ_F;
               invf = DZ_INVF;
               break;
            case VR_INDEX:
               f = VR_F;
               invf = VR_INVF;
               break;
            case SW_INDEX:
               f = VR_F;
               invf = VR_INVF;
               break;
            case CZ_INDEX:
               f = DZ_F;
               invf = DZ_INVF;
               break;
            case DR_INDEX:
               f = DR_F;
               invf = DR_INVF;
               break;
            case PH_INDEX:
            case NP_INDEX:
               f = PH_F;
               invf = PH_INVF;
               break;
            case KD_INDEX:
               f = KD_F;
               invf = KD_INVF;
               break;
            case RH_INDEX:
               f = RH_F;
               invf = RH_INVF;
            break;
            }
         volume[n_volume]->h.f = f;
         volume[n_volume]->h.invf = invf;
         volume[n_volume]->sweep[n_sweep]->h.f = f;
         volume[n_volume]->sweep[n_sweep]->h.invf = invf;
//         strcpy(volume[n_volume]->h.type_str, var_name);
         
         for (n_ray = 0; n_ray < NRAYS; n_ray++)
            {
            ray = RSL_new_ray(NBINS);
            date_tmp[0] = date_str[0];
            date_tmp[1] = date_str[1];
            date_tmp[2] = date_str[2];
            date_tmp[3] = date_str[3];
            date_tmp[4] = '\0';
            ray->h.year = atoi(date_tmp);
            date_tmp[0] = date_str[4];
            date_tmp[1] = date_str[5];
            date_tmp[2] = '\0';
            ray->h.month = atoi(date_tmp);
            date_tmp[0] = date_str[5];
            date_tmp[1] = date_str[7];
            date_tmp[2] = '\0';
            ray->h.day = atoi(date_tmp);
            date_tmp[0] = time_str[0];
            date_tmp[1] = time_str[1];
            date_tmp[2] = '\0';
            ray->h.hour = atoi(date_tmp);
            date_tmp[0] = time_str[2];
            date_tmp[1] = time_str[3];
            date_tmp[2] = '\0';
            ray->h.minute = atoi(date_tmp);
            date_tmp[0] = time_str[4];
            date_tmp[1] = time_str[5];
            date_tmp[2] = '\0';
            ray->h.sec = atoi(date_tmp);
                        ray->h.ray_num = n_ray + 1;
                        
            ray->h.gate_size = r_step;
            ray->h.range_bin1 = r_start;
            ray->h.elev_num = volume[n_volume]->sweep[n_sweep]->h.sweep_num;
	    ray->h.elev = elev;
            ray->h.azimuth = (azAngleIni[n_ray] + azAngleEnd[n_ray])/2;
            ray->h.prf = prf;
            ray->h.prf2 = prf2;
            ray->h.pulse_width = pulsewidth;
            ray->h.wavelength = wavelength;
            ray->h.beam_width = (verti_beam + horiz_beam)/2;
            ray->h.nyq_vel = nyq_vel;
            ray->h.f = f;
            ray->h.invf = invf;
            ray->h.nbins = NBINS;
            if (H5T_FLOAT == class_id)
               {
               for (n_bin = 0; n_bin < NBINS; n_bin++)
                  {
                  vtemp = image32[n_bin + NBINS*n_ray]*gain + offset;
                  if ((vtemp == nodata) || (vtemp == undetect) || (vtemp == offset))
                     {
                     vtemp = NOECHO;
                     }
                  ray->range[n_bin] = invf(vtemp);
                  }
               }
            else
               {
               for (n_bin = 0; n_bin < NBINS; n_bin++)
                  {
                  vtemp = image[n_bin + NBINS*n_ray]*gain + offset;
                  if ((vtemp == nodata) || (vtemp == undetect) || (vtemp == offset))
                     {
                     vtemp = NOECHO;
                     }
                  ray->range[n_bin] = invf(vtemp);
                  }
               }
            volume[n_volume]->sweep[n_sweep]->ray[n_ray] = ray;
            }
         
//         printf("volume %s\t MAX = %10.4f   MIN = %10.4f  GAIN = %10.4f  OFFSET = %10.4f\n", var_name, vmax, vmin, gain, offset);
         if (H5T_INTEGER == class_id)
            {   
            free(image);
            }
         else
            {
            free(image32);
            }
         H5Gclose(vol);
         }
      H5Gclose(scan);
      }

   for (n_volume = 0; n_volume < MAX_RADAR_VOLUMES; n_volume++)
      {
      if (-1 != use_volume[n_volume])
         {
         radar->v[use_volume[n_volume]] = volume[n_volume];
         valid_vol++;
         }
      }
   /*
   for (n_volume = valid_vol; n_volume < MAX_RADAR_VOLUMES; n_volume++)
      {
      radar->v[n_volume++] = NULL;
      }
   */
   radar->h.nvolumes = valid_vol;
   
   return radar;
  }

int return_type_var(char *varname)
   {
   if (0 == strncmp(varname, "TH", strlen(varname)))
      {
      return DZ_INDEX;
      }
   if (0 == strncmp(varname, "DBZH", strlen(varname)))
      {
      return CZ_INDEX;
      }
   if (0 == strncmp(varname, "VRAD", strlen(varname)))
      {
      return VR_INDEX;
      }
   if (0 == strncmp(varname, "WRAD", strlen(varname)))
      {
      return SW_INDEX;
      }
   return -1;
   }


Radar *RSL_hdf5_to_radar(char *infile)
  {
  hid_t file, group, vol, attr, dataspace, memspace, how,
    header, elev0, elev1, azim0, azim1, what, memtype, where,
    time;
  int status;
  int variavel;
  
  int n_sweep, n_ray, n_bin, n_volume;
  char grp_name[MAX_HDF_STR];
  int rank;
  /*
  int volindx[DEFINED_VOLUMES] = {CZ_INDEX, DZ_INDEX, VR_INDEX, 
				  SW_INDEX, DR_INDEX, PH_INDEX, 
				  NP_INDEX, KD_INDEX, RH_INDEX};
  */

  int volindx[DEFINED_VOLUMES];

  int vartype = -1;
  
  
  char moments[MAX_RADAR_VOLUMES][MAX_HDF_STR];
  char date_tmp[MAX_HDF_STR];
  
  int pulse, prf;
  double r_step, r_sample, r_start, horiz_beam, verti_beam,
    height, lat, lon, min, sec, wavelength;
  float r_min, r_max, scale_factor, elev_header;
  char arquivo[MAX_HDF_STR];
  
  hsize_t *dims = NULL;
  int *data_out = NULL;
  double *ray_elev0 = NULL, *ray_elev1 = NULL;
  double *ray_azim0 = NULL, *ray_azim1 = NULL;
   unsigned long long *timestamp = NULL;
   char *date_time = NULL, *var_name = NULL,
      *template_name = NULL, *data_size = NULL;
   Radar *radar = NULL;
   Sweep *sweep = NULL;
   Volume *volume = NULL;
   Ray *ray = NULL;
   struct tm *cal_time = NULL;

   static float (*f)(Range x);
   static Range (*invf)(float x);

   int NSWEEPS, NVOLUMES, NRAYS, NBINS;

   time_t scan_time;

   char version_str[MAX_HDF_STR];

#ifndef VERSION1_8_0
//   herr_t (*error_stack)(void*);
   
   
   /* Turn off error handling */
//   H5Eset_auto(error_stack, NULL);
#endif
   
   strcpy(arquivo, infile);
   
   file = H5Fopen(arquivo, H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0)
      {
      printf("ERRO\n");
      return (-1);
      }


   for (n_volume=0; n_volume<MAX_RADAR_VOLUMES; n_volume++)
      {
      sprintf(moments[n_volume], "moment_%d\0", n_volume);
      }
   
   radar = RSL_new_radar(MAX_RADAR_VOLUMES);

   /*Read information about the beam width*/
#ifdef VERSION1_8_0   
   how = H5Gopen(file, "/how", H5P_DEFAULT);
   what = H5Gopen(file, "/what", H5P_DEFAULT);
   where = H5Gopen(file, "/where", H5P_DEFAULT);
#else
   how = H5Gopen(file, "/how");
   what = H5Gopen(file, "/what");
   where = H5Gopen(file, "/where");
#endif


   /*verifica se e HDF ODIM ou GAMIC*/
   attr = H5Aopen(what, "version", H5P_DEFAULT);
   memtype = H5Aget_type(attr);
   memset(version_str, 0, sizeof(version_str));
   status = H5Aread(attr, memtype, &version_str);
   H5Tclose(memtype);
   H5Aclose(attr);
   if (0 == strncmp(version_str, "H5Rad 2.1", 9))
      {
      H5Gclose(how);
      H5Gclose(what);
      H5Gclose(where);
      H5Fclose(file);
      return RSL_hdf5_ODIM_to_radar(arquivo);
      }
   
   
   attr = H5Aopen(how, "azimuth_beam", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &horiz_beam);
   H5Aclose(attr);

   attr = H5Aopen(how, "elevation_beam", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &verti_beam);
   H5Aclose(attr);

   attr = H5Aopen(where, "height", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &height);
   H5Aclose(attr);

   attr = H5Aopen(where, "lat", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &lat);
   H5Aclose(attr);

   attr = H5Aopen(where, "lon", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_DOUBLE, &lon);
   H5Aclose(attr);

   /*reading date/time information*/
   attr = H5Aopen(what, "date", H5P_DEFAULT);
   memtype = H5Tcopy(H5T_C_S1);
   /*string size in BYTES*/
   H5Tset_size(memtype, H5T_VARIABLE);
   status = H5Aread(attr, memtype, &date_time);
   H5Aclose(attr);
   /*reading template name*/

#if 0
   attr = H5Aopen(how, "template_name", H5P_DEFAULT);
   memtype = H5Tcopy(H5T_C_S1);
   //   H5Tset_size(memtype, MAX_HDF_STR);
   H5Tset_size(memtype, H5T_VARIABLE);

   status = H5Aread(attr, memtype, &template_name);
   H5Aclose(attr);
#endif
   H5Tclose(memtype);
   H5Gclose(how);
   H5Gclose(what);
   H5Gclose(where);
 
   date_tmp[0] = date_time[0];
   date_tmp[1] = date_time[1];
   date_tmp[2] = date_time[2];
   date_tmp[3] = date_time[3];
   date_tmp[4] = '\0';
   radar->h.year = atoi(date_tmp);
   date_tmp[0] = date_time[5];
   date_tmp[1] = date_time[6];
   date_tmp[2] = '\0';
   radar->h.month = atoi(date_tmp);
   date_tmp[0] = date_time[8];
   date_tmp[1] = date_time[9];
   date_tmp[2] = '\0';
   radar->h.day = atoi(date_tmp);
   date_tmp[0] = date_time[11];
   date_tmp[1] = date_time[12];
   date_tmp[2] = '\0';
   radar->h.hour = atoi(date_tmp);
   date_tmp[0] = date_time[14];
   date_tmp[1] = date_time[15];
   date_tmp[2] = '\0';
   radar->h.minute = atoi(date_tmp);
   date_tmp[0] = date_time[17];
   date_tmp[1] = date_time[18];
   date_tmp[2] = '\0';
   radar->h.sec = atoi(date_tmp);
   
   free(date_time);
   free(template_name);

   /*writes lat/lon and height info to the radar header*/
   radar->h.height = height;

   radar->h.latd = (int) lat;
   radar->h.lond = (int) lon;
   lat = fabs(lat);
   lon = fabs(lon);
   radar->h.latm = (int) ((lat - abs(radar->h.latd)) * 60);
   radar->h.lonm = (int) ((lon - abs(radar->h.lond)) * 60);
   radar->h.lats = (int) ((((lat - abs(radar->h.latd)) * 60) -
                           radar->h.latm) * 60);
   radar->h.lons = (int) ((((lon - abs(radar->h.lond)) * 60) -
                           radar->h.lonm) * 60);
   /*
   - locates the number of scans
   - locates the first scan and finds the total volume number
   */
#ifdef VERSION1_8_0   
   what = H5Gopen(file, "/what", H5P_DEFAULT);
#else
   what = H5Gopen(file, "/what");
#endif
   attr = H5Aopen(what, "sets", H5P_DEFAULT);
   status = H5Aread(attr, H5T_NATIVE_INT, &NSWEEPS);
   H5Aclose(attr);
   H5Gclose(what);
   
   for (n_sweep=0; n_sweep < NSWEEPS; n_sweep++)
      {
      memset(grp_name, 0, sizeof(grp_name));
      sprintf(grp_name, "scan%d", n_sweep);
#ifdef VERSION1_8_0
      group = H5Gopen(file, grp_name, H5P_DEFAULT);
#else
      group = H5Gopen(file, grp_name);
#endif
      if (group >= 0)
         {
         NVOLUMES = 0;
         for (n_volume = 0; n_volume < DEFINED_VOLUMES; n_volume++)
            {
#ifdef VERSION1_8_0
            vol = H5Dopen(group, moments[n_volume], H5P_DEFAULT);
#else
            vol = H5Dopen(group, moments[n_volume]);
#endif
            if (vol >= 0)
               {
               attr = H5Aopen(vol, "moment", H5P_DEFAULT);
               memtype = H5Tcopy(H5T_C_S1);
	       //               H5Tset_size(memtype, MAX_HDF_STR);
	       H5Tset_size(memtype, H5T_VARIABLE);

               status = H5Aread(attr, memtype, &var_name);
               H5Aclose(attr);
               vartype = return_vartype(var_name);
               
               if (vartype > -1)
                  {
                  volindx[NVOLUMES] = vartype;
                  NVOLUMES++;
                  }
               free(var_name);
               }
            }
         H5Gclose(group);
         break;
         }
      }
   
   for (n_volume = 0; n_volume < NVOLUMES; n_volume++)
      {
      volume = RSL_new_volume(NSWEEPS);
      switch (volindx[n_volume])
         {
         case DZ_INDEX:
            f = DZ_F;
            invf = DZ_INVF;
            break;
         case VR_INDEX:
            f = VR_F;
            invf = VR_INVF;
            break;
         case SW_INDEX:
            f = VR_F;
            invf = VR_INVF;
            break;
         case CZ_INDEX:
            f = DZ_F;
            invf = DZ_INVF;
            break;
         case DR_INDEX:
            f = DR_F;
            invf = DR_INVF;
            break;
         case PH_INDEX:
         case NP_INDEX:
            f = PH_F;
            invf = PH_INVF;
            break;
         case KD_INDEX:
            f = KD_F;
            invf = KD_INVF;
            break;
         case RH_INDEX:
            f = RH_F;
            invf = RH_INVF;
            break;
         }
      volume->h.f = f;
      volume->h.invf = invf;
      for (n_sweep=0; n_sweep < NSWEEPS; n_sweep++)
         {
         memset(grp_name, 0, sizeof(grp_name));
         sprintf(grp_name, "scan%d", n_sweep);
#ifdef VERSION1_8_0         
         group = H5Gopen(file, grp_name, H5P_DEFAULT);
         if (group < 0)
            continue;
         
         vol = H5Dopen(group, moments[n_volume], H5P_DEFAULT);
         if (vol < 0)
            continue;

         header = H5Dopen(group, "ray_header", H5P_DEFAULT);
         if (header < 0)
            {
            H5Dclose(vol);
            continue;
            }

         how = H5Gopen(group, "how", H5P_DEFAULT);
         if (how < 0)
            {
            H5Dclose(vol);
            H5Dclose(header);
            continue;
            }
#else
         group = H5Gopen(file, grp_name);
         if (group < 0)
            continue;
         
         vol = H5Dopen(group, moments[n_volume]);
         if (vol < 0)
            continue;
         
         header = H5Dopen(group, "ray_header");
         if (header < 0)
            {
            H5Dclose(vol);
            continue;
            }
         
         how = H5Gopen(group, "how");
         if (how < 0)
            {
            H5Dclose(vol);
            H5Dclose(header);
            continue;
            }
#endif            
         /*how to read the data...*/
         dataspace = H5Dget_space (vol); /* dataspace handle */
         rank = H5Sget_simple_extent_ndims (dataspace);
         dims = (hsize_t *) calloc(rank, sizeof(hsize_t));
         H5Sget_simple_extent_dims (dataspace, dims, NULL);
         NRAYS = dims[0];
         NBINS = dims[1];
         
         /*output buffer*/
         data_out = (int *) calloc(NRAYS*NBINS, sizeof(int));
         status = H5Dread(vol, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, data_out);
         
         
         /*how to read the ray header...*/
         dataspace = H5Dget_space (header);
         rank = H5Sget_simple_extent_ndims (dataspace);
         dims = (hsize_t *) calloc(rank, sizeof(hsize_t));  
         H5Sget_simple_extent_dims (dataspace, dims, NULL);
         /*output buffer*/
         
         ray_elev0 = (double *) calloc(NRAYS, sizeof(double));
         ray_elev1 = (double *) calloc(NRAYS, sizeof(double));
         ray_azim0 = (double *) calloc(NRAYS, sizeof(double));
         ray_azim1 = (double *) calloc(NRAYS, sizeof(double));
         timestamp = (unsigned long long *) calloc(NRAYS, sizeof(unsigned long long));
         
         azim0 = H5Tcreate(H5T_COMPOUND, sizeof(double));
         status = H5Tinsert(azim0, "azimuth_start", 0,
                            H5T_NATIVE_DOUBLE);
         status = H5Dread(header, azim0, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, ray_azim0);
         H5Tclose(azim0);
         
         
         elev0 = H5Tcreate(H5T_COMPOUND, sizeof(double));
         status = H5Tinsert(elev0, "elevation_start", 0,
                            H5T_NATIVE_DOUBLE);
         status = H5Dread(header, elev0, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, ray_elev0);
         H5Tclose(elev0);
         
         elev1 = H5Tcreate(H5T_COMPOUND, sizeof(double));
         status = H5Tinsert(elev1, "elevation_stop", 0,
                            H5T_NATIVE_DOUBLE);
         status = H5Dread(header, elev1, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, ray_elev1);
         H5Tclose(elev1);

         time = H5Tcreate(H5T_COMPOUND, sizeof(unsigned long long));
         status = H5Tinsert(time, "timestamp", 0,
                            H5T_NATIVE_ULLONG);
         status = H5Dread(header, time, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, timestamp);
         H5Tclose(time);

         
         /*reading sweep and moment info*/
         attr = H5Aopen(how, "PRF", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_INT, &prf);
         H5Aclose(attr);


	 if (H5Aexists(how, "pulse_width_mks"))
	   {
	   attr = H5Aopen(how, "pulse_width_mks", H5P_DEFAULT);
	   }
	 else
	   {
	   attr = H5Aopen(how, "pulse_width_us", H5P_DEFAULT);
	   }

         status = H5Aread(attr, H5T_NATIVE_INT, &pulse);
         H5Aclose(attr);

         attr = H5Aopen(how, "radar_wave_length", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_DOUBLE, &wavelength);
         H5Aclose(attr);

         attr = H5Aopen(how, "range_samples", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_DOUBLE, &r_sample);
         H5Aclose(attr);

         attr = H5Aopen(how, "range_start", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_DOUBLE, &r_start);
         H5Aclose(attr);

         attr = H5Aopen(how, "range_step", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_DOUBLE, &r_step);
         H5Aclose(attr);

         attr = H5Aopen(vol, "dyn_range_min", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_FLOAT, &r_min);
         H5Aclose(attr);

         attr = H5Aopen(vol, "dyn_range_max", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_FLOAT, &r_max);
         H5Aclose(attr);

         attr = H5Aopen(how, "elevation", H5P_DEFAULT);
         status = H5Aread(attr, H5T_NATIVE_FLOAT, &elev_header);
         H5Aclose(attr);

         attr = H5Aopen(vol, "format", H5P_DEFAULT);
         memtype = H5Tcopy(H5T_C_S1);
	 //         H5Tset_size(memtype, MAX_HDF_STR);
	 H5Tset_size(memtype, H5T_VARIABLE);

         status = H5Aread(attr, memtype, &data_size);
         H5Aclose(attr);

//         printf("RANGE = %d %f %f\n", n_volume, r_min, r_max);
         
         
	 // r_max = 50;
	 //r_min = -50;


         if (strncmp(data_size, "UV16", 4))
            {   
            scale_factor = (r_max - r_min) / 254;
            }
         else
            {
            scale_factor = (r_max - r_min) / 65534;
            }
         
         
         volume->sweep[n_sweep] = RSL_new_sweep(NRAYS);
         volume->sweep[n_sweep]->h.sweep_num = n_sweep + 1;
         volume->sweep[n_sweep]->h.beam_width = (verti_beam + horiz_beam)/2;
         volume->sweep[n_sweep]->h.vert_half_bw = verti_beam / 2;
         volume->sweep[n_sweep]->h.horz_half_bw = horiz_beam / 2;
	 //         printf("ELEV: %f\n", elev_header);
         for (n_ray = 0; n_ray < NRAYS; n_ray++)
            {
            ray = RSL_new_ray(NBINS);
            /* timestamp info in unix_time * 1000000 */
            scan_time = (time_t) (timestamp[n_ray] / 1000000);
            cal_time = gmtime(&scan_time);
            ray->h.year = cal_time->tm_year + 1900;
            ray->h.month = cal_time->tm_mon + 1;
            ray->h.day = cal_time->tm_mday;
            ray->h.hour = cal_time->tm_hour;
            ray->h.minute = cal_time->tm_min;
            ray->h.sec = cal_time->tm_sec;
            
            ray->h.ray_num = n_ray + 1;
            ray->h.gate_size = r_step;
	    ray->h.gate_size = r_step*r_sample;
            ray->h.range_bin1 = r_start;
            ray->h.elev_num = volume->sweep[n_sweep]->h.sweep_num;
	    ray->h.elev = (float) ((ray_elev0[n_ray] + ray_elev1[n_ray])/2);
	    //ray->h.elev = elev_header;
            ray->h.azimuth = ray_azim0[n_ray];
            ray->h.prf = prf;
            ray->h.pulse_width = pulse;
            ray->h.wavelength = wavelength;
            ray->h.beam_width = (verti_beam + horiz_beam)/2;
            /*
            ray->h.nyq_vel = Scount[i].uvv;
            */
            ray->h.f = f;
            ray->h.invf = invf;
            ray->h.nbins = NBINS;
            ray->range[NBINS - 1] = invf(NOECHO);
            for (n_bin = 0; n_bin < NBINS - 1; n_bin++)
               {
		 //printf("%d %d %d\n", n_ray, n_bin, data_out[n_ray*NBINS + n_bin]);
               if (data_out[n_ray*NBINS + n_bin] == 0)
                  {
                  ray->range[n_bin] = invf(NOECHO);
                  continue;
                  }
               
               ray->range[n_bin] = invf(((((float) data_out[n_ray*NBINS + n_bin]) - 1) * scale_factor + r_min));
               
#if 0
               switch (volindx[n_volume])
                  {
                  case DZ_INDEX:
                  case CZ_INDEX:
                     ray->range[n_bin] = invf(((float) data_out[n_ray*NBINS + n_bin])/2.0 - 32);
                     break;
                  case VR_INDEX:
                     /*
                     range is 1 - 255, zero is zero...
                     1 == r_min and 255 == r_max
                     0 == 0
                     */
                     ray->range[n_bin] = invf((((float) data_out[n_ray*NBINS + n_bin]) * scale_factor + r_min));      
                     break;
                  case SW_INDEX:
                     
                     ray->range[n_bin] = invf(((float) data_out[n_ray*NBINS + n_bin])/255.0);
                     break;
                  }
#endif            
               }
            volume->sweep[n_sweep]->ray[n_ray] = ray;
            }
         volume->sweep[n_sweep]->h.nrays = NRAYS;
	 volume->sweep[n_sweep]->h.elev = elev_header;
         volume->sweep[n_sweep]->h.f = f;
         volume->sweep[n_sweep]->h.invf = invf;
         
         /*Close dataspace*/
         H5Sclose(dataspace);
         free(dims);
         free(data_out);
         free(ray_elev0);
         free(ray_elev1);
         free(ray_azim0);
         free(ray_azim1);
         free(timestamp);
         H5Dclose(vol);
         H5Dclose(header);
         H5Gclose(how);
         }
      radar->v[volindx[n_volume]] = volume;
      }
   H5Fclose(file);

   return radar;
   }
#else
Radar *RSL_hdf5_to_radar(char *infile)
{
  return NULL;
}
#endif

#ifdef HAVE_LIBTSDISTK
/******************************************************************
   Reads one volume scan from a HDF file into a RSL radar structure.

  -----------------------------------------------------------------
	 Libraries required for execution of this code :
      -ltsdistk                    : tsdis toolkit
      -lmfhdf -ldf -ljpeg -lz      : HDF
      -lrsl                        : rsl
      -lm                          : C math

  -----------------------------------------------------------------
*******************************************************************/

/* TSDIS toolkit function and structure definitions. */
#include "IO.h"
#include "IO_GV.h"
/* Parameter definitions for 1B-51 and 1C-51 HDF
	 file handling applications using the TSDIS toolkit. */
#include "toolkit_1BC-51_appl.h"
#define MISSING_VAL 0



/*************************************************************/
/*                                                           */
/*                Function Prototypes                        */
/*                                                           */
/*************************************************************/
void RayFillFrom1B51(Ray *ray, int16 *rayData, PARAMETER_DESCRIPTOR *parmDesc);
void RayFillFrom1C51(Ray *ray, int vindex, int16 *rayData, int8 *rayMaskData,
										 PARAMETER_DESCRIPTOR *parmDesc, float calibr);
static void Ray_headerFill(Ray *ray, L1B_1C_GV *gvl1, VosSize *vs, 
													 int pindex, int tk_sindex, int rindex);
Ray *RayBuild(L1B_1C_GV *gvl1, VosSize *vs, float calibr, 
							int vindex, int pindex, int sindex, int rindex);
static void Sweep_headerFill(Sweep *sweep, SENSORS *sensor, int sindex, int nrays);
Sweep *SweepBuild(L1B_1C_GV *gvl1, VosSize *vs, float calibr, 
									int vindex, int pindex, int sindex);
static void Volume_headerFill(Volume *volume, char *parmDesc, int vindex,
											 int nsweeps, float calibr);
Volume *VolumeBuild(L1B_1C_GV *gvl1, VosSize *vs, float calibr, 
										int vindex, int pindex);
int parmIdentify(char *parmName);
static void Radar_headerFill(Radar *radar, L1B_1C_GV *gvl1);
Radar *RadarBuild(L1B_1C_GV *gvl1, VosSize *vs, float zCal);
int commentsRead(VosSize *vs, float *zCal, char *comments, int productID);
static L1B_1C_GV *GVL1Build(IO_HANDLE *granuleHandle, int vosNum,
														VosSize *vs);
int metaDataRead(Radar *radar, IO_HANDLE *granuleHandle);
static int hdfFileOpen(char *infile, IO_HANDLE *granuleHandle, 
											 char *hdfFileName, int *vosNum);
Radar *RSL_hdf_to_radar(char *infile);


/* Toolkit memory management functions. */
extern void TKfreeGVL1(L1B_1C_GV *gvl1);
extern int8 ***TKnewParmData1byte(int nsweep, int nray, int ncell);
extern int16 ***TKnewParmData2byte(int nsweep, int nray, int ncell);
extern PARAMETER *TKnewGVL1parm(void);
extern L1B_1C_GV *TKnewGVL1(void);


static float (*f)(Range x);
static Range (*invf)(float x);

extern int radar_verbose_flag;


static void ymd(int jday, int yy, int *mm, int *dd);

static int daytab[2][13] = {
  {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365},
  {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}
};

static void ymd(int jday, int year, int *mm, int *dd)
{
  /*  Input: jday, yyyy */
  /* Output: mm, dd */
  int leap;
  int i;

  leap = (year%4 == 0 && year%100 != 0) || year%400 == 0;
  for (i=0; daytab[leap][i]<jday; i++) continue;
  *mm = i;
  i--;
  *dd = jday - daytab[leap][i];
}

/*************************************************************/
/*                                                           */
/*                       RayFillFrom1B51                     */
/*                                                           */
/*************************************************************/
void RayFillFrom1B51(Ray *ray, int16 *rayData, PARAMETER_DESCRIPTOR *parmDesc)
{
/* 
	 Fill the RSL bin slots of one ray of any volume using the corresponding
	 ray data from a 1B-51 HDF file.
*/
	int j; 
	
	for (j=0; j<ray->h.nbins; j++)
	{
		if (rayData[j] <= AP_VALUE)  /* Handle anomalous condition flags. */
		{
			if (rayData[j] == NO_VALUE) ray->range[j] = invf((float)BADVAL);
			else if (rayData[j] == RNG_AMBIG_VALUE) ray->range[j] = invf((float)RFVAL);
			else if (rayData[j] == NOECHO_VALUE) ray->range[j] = invf((float)NOECHO);
			else ray->range[j] = invf((float)APFLAG);
		}
		else  /* Valid data value */
		{
			ray->range[j] = invf( (rayData[j] - parmDesc->offsetFactor) / 
				                     parmDesc->scaleFactor );
		}
	} /* end for (j=0... */
}

/*************************************************************/
/*                                                           */
/*                         RayFillFrom1C51                   */
/*                                                           */
/*************************************************************/
void RayFillFrom1C51(Ray *ray, int vindex, int16 *rayData, int8 *rayMaskData,
										 PARAMETER_DESCRIPTOR *parmDesc, float calibr)
{
/* 
	 Fill the RSL bin slots of one ray of a CZ or a DZ volume using the 
	 corresponding ray data and ray_mask data from a 1C-51 HDF file.
*/
	int j; 

	for (j=0; j<ray->h.nbins; j++)
	{
		if (rayData[j] <= AP_VALUE)  /* Handle anomalous condition flags. */
		{
			if (rayData[j] == NO_VALUE) ray->range[j] = invf((float)BADVAL);
			else if (rayData[j] == RNG_AMBIG_VALUE) ray->range[j] = invf((float)RFVAL);
			else if (rayData[j] == NOECHO_VALUE) ray->range[j] = invf((float)NOECHO);
			else ray->range[j] = invf((float)APFLAG);
		}
		else  /* Valid data value */
		{
			if ((vindex == CZ_INDEX) || (vindex == CD_INDEX))
			{
				if (rayMaskData[j] == 1)
				  ray->range[j] = invf((float)BADVAL);
				else
				  ray->range[j] = invf( (rayData[j] - parmDesc->offsetFactor) / 
					 parmDesc->scaleFactor );
			}  /* end if (vindex == CZ_INDEX) */
			else if (vindex == DZ_INDEX)
			{
				ray->range[j] = invf( 
							((rayData[j] - parmDesc->offsetFactor) / parmDesc->scaleFactor) -
								calibr + X * rayMaskData[j] );
				
			}  /* end else if DZ_INDEX */
			else if (vindex == ZD_INDEX)
			{
				ray->range[j] = invf( 
						((rayData[j] - parmDesc->offsetFactor) / parmDesc->scaleFactor) + 
							X * rayMaskData[j] );
			}  /* end else if ZD_INDEX */
			else
			  fprintf(stderr, "RayFillFrom1C51(): Illegal volume index..\n");
		} /* else valid data value */
	}  /* end for (j=0; ... */
}

/*************************************************************/
/*                                                           */
/*                      Ray_headerFill                       */
/*                                                           */
/*************************************************************/
void Ray_headerFill(Ray *ray, L1B_1C_GV *gvl1, VosSize *vs, 
										int pindex, int tk_sindex, int rindex)
{
	ray->h.year = (int)gvl1->volDes.year;
	/* Get calendar date (month, day) from (year, Julian day) */
	ymd((int)gvl1->sensor.rayInfoInteger[tk_sindex][rindex][1],
		ray->h.year, &ray->h.month, &ray->h.day);
	ray->h.hour = (int)gvl1->sensor.rayInfoInteger[tk_sindex][rindex][2];
	ray->h.minute = (int)gvl1->sensor.rayInfoInteger[tk_sindex][rindex][3];
	ray->h.sec = (float)(gvl1->sensor.rayInfoInteger[tk_sindex][rindex][4] + 
											 gvl1->sensor.rayInfoInteger[tk_sindex][rindex][5]/1000.0);

	ray->h.azimuth = gvl1->sensor.rayInfoFloat[tk_sindex][rindex][0];
	ray->h.ray_num = rindex + 1;
	ray->h.elev = gvl1->sensor.rayInfoFloat[tk_sindex][rindex][1]; /* degrees */
	ray->h.elev_num = tk_sindex + 1;
	ray->h.gate_size = (int)
	(gvl1->sensor.parm[pindex]->cellRangeVector.distanceToCell[2] -
	 gvl1->sensor.parm[pindex]->cellRangeVector.distanceToCell[1]); /*meters*/

	ray->h.range_bin1 = (int)
	          (gvl1->sensor.parm[pindex]->cellRangeVector.distanceToCell[0] - 
	                                        0.5 * ray->h.gate_size); /* meters */

	ray->h.vel_res = MISSING_VAL;           /* ?? */
	/* Sweeps/min */
	ray->h.sweep_rate = (float)(gvl1->sensor.radarDesc.nomScanRate / 6.0);  
	ray->h.prf = (int)gvl1->sensor.rayInfoFloat[tk_sindex][rindex][3];
	ray->h.azim_rate = (float)gvl1->sensor.radarDesc.nomScanRate;
	ray->h.fix_angle = (float)gvl1->sensor.sweepInfo[tk_sindex].fixedAngle;
	ray->h.pulse_count = (int)gvl1->sensor.rayInfoFloat[tk_sindex][rindex][2];
	/* Pulse width (microsec) */
	ray->h.pulse_width = (float)(gvl1->sensor.parm[pindex]->parmDesc.pulseWidth /
															 300.0);
	ray->h.beam_width = (float)gvl1->sensor.radarDesc.horBeamWidth;
	/* Carrier freq (GHz) */
	ray->h.frequency = (float)gvl1->sensor.radarDesc.frequency1;
	/* wavelength (m) */
	if (ray->h.frequency != 0.0)
	  ray->h.wavelength = (RSL_SPEED_OF_LIGHT / ray->h.frequency) * 1.0e-9;
	else
	  ray->h.wavelength = 0.0;
	ray->h.nyq_vel = (float) (ray->h.prf * ray->h.wavelength / 4.0);
	if (ray->h.prf != 0)
	  ray->h.unam_rng = (float) RSL_SPEED_OF_LIGHT / (2.0 * ray->h.prf * 1000.0);
	else
	  ray->h.unam_rng = (float) 0.0;
	ray->h.nbins = vs->tk.ncell[tk_sindex][pindex];
	ray->h.f = f;
	ray->h.invf = invf;
}

/*************************************************************/
/*                                                           */
/*                          RayBuild                         */
/*                                                           */
/*************************************************************/
Ray *RayBuild(L1B_1C_GV *gvl1, VosSize *vs, float calibr, 
							int vindex, int pindex, int tk_sindex, int rindex)
{
	Ray *ray;
	
	/* Create a Ray structure. */
	ray = RSL_new_ray(vs->tk.ncell[tk_sindex][pindex]);
	if (ray == NULL)
	{
		perror("RayBuild(): RSL_new_ray failed\n");
		return(NULL);
	}
	Ray_headerFill(ray, gvl1, vs, pindex, tk_sindex, rindex);
	/* Is this a 1C-51 file? */
	if ((strcmp(gvl1->sensor.parm[pindex]->parmDesc.parmName, "QCZ") == 0) ||
			(strcmp(gvl1->sensor.parm[pindex]->parmDesc.parmName, "QCZDR") == 0))
	  RayFillFrom1C51(ray, vindex,
						gvl1->sensor.parm[pindex]->parmData2byte[tk_sindex][rindex],
						gvl1->sensor.parm[pindex-1]->parmData1byte[tk_sindex][rindex],
						&gvl1->sensor.parm[pindex]->parmDesc, calibr);
	else /* 1B-51 file */
	  RayFillFrom1B51(ray,
						gvl1->sensor.parm[pindex]->parmData2byte[tk_sindex][rindex],
						&gvl1->sensor.parm[pindex]->parmDesc);
	return(ray);
}

/*************************************************************/
/*                                                           */
/*                     Sweep_headerFill                      */
/*                                                           */
/*************************************************************/
void Sweep_headerFill(Sweep *sweep, SENSORS *sensor, int tk_sindex, int nrays)
{
/*	sweep->h.sweep_num filled in VolumeBuild() */
	sweep->h.elev = sensor->sweepInfo[tk_sindex].fixedAngle;
	sweep->h.beam_width = sensor->radarDesc.horBeamWidth;
	sweep->h.horz_half_bw = sensor->radarDesc.horBeamWidth / 2.0;
	sweep->h.vert_half_bw = sensor->radarDesc.verBeamWidth / 2.0;
	sweep->h.nrays = sensor->sweepInfo[tk_sindex].numRays;
	sweep->h.f = f;
	sweep->h.invf = invf;
}

/*************************************************************/
/*                                                           */
/*                        SweepBuild                         */
/*                                                           */
/*************************************************************/
Sweep *SweepBuild(L1B_1C_GV *gvl1, VosSize *vs, float calibr, 
									int vindex, int pindex, int tk_sindex)
{
	int rindex;
	Sweep *sweep;
	
	/* Create a Sweep structure. */
	sweep = RSL_new_sweep(vs->rsl.maxNray);
	if (sweep == NULL)
	{
		perror("SweepBuild(): RSL_new_sweep failed\n");
		return(NULL);
	}
	/* Initialize the Sweep_header values. */
	Sweep_headerFill(sweep, &gvl1->sensor, tk_sindex, vs->rsl.maxNray);

	/* Loop to fill each of the rays of this rsl sweep structure. */
	for (rindex=0; rindex<vs->tk.nray[tk_sindex]; rindex++)
		sweep->ray[rindex] = RayBuild(gvl1, vs, calibr, vindex, pindex,
																	tk_sindex, rindex);
	return(sweep);
}

/*************************************************************/
/*                                                           */
/*                    Volume_headerFill                      */
/*                                                           */
/*************************************************************/
void Volume_headerFill(Volume *volume, char *parmDesc, int vindex,
											 int nsweeps, float calibr)
{
	if (vindex == DZ_INDEX)
	  volume->h.type_str = strdup("Reflectivity");
	else if (vindex == ZD_INDEX)
	  volume->h.type_str = strdup("Differential Reflectivity");
	else volume->h.type_str = strdup(parmDesc);
	volume->h.f = f;
	volume->h.invf = invf;
	volume->h.nsweeps = nsweeps;
	volume->h.calibr_const = calibr;
}

/*************************************************************/
/*                                                           */
/*                        VolumeBuild                        */
/*                                                           */
/*************************************************************/
Volume *VolumeBuild(L1B_1C_GV *gvl1, VosSize *vs, float calibr, 
										int vindex, int pindex)
{		
	Volume *v;
	int sindex, tk_sindex;
	extern int *rsl_qsweep; /* See RSL_read_these_sweeps in volume.c */
	extern int rsl_qsweep_max;
	
	/* Create a Volume structure. */
	v = RSL_new_volume(vs->tk.nsweep);
	if (v == NULL)
	{
		perror("VolumeBuild(): RSL_new_volume failed\n");
		return(NULL);
	}
	/* Initialize the Volume_header values. */
	Volume_headerFill(v, gvl1->sensor.parm[pindex]->parmDesc.parmDesc, 
										vindex, vs->tk.nsweep, calibr);
	if (radar_verbose_flag)
	  fprintf(stderr, "RSL volume type: %s\n", v->h.type_str);

	/* Build each of the sweeps of this radar volume structure. */
	sindex = -1;
	for (tk_sindex=0; tk_sindex<vs->tk.nsweep; tk_sindex++)
	{
	  if (rsl_qsweep != NULL) {
		if (tk_sindex > rsl_qsweep_max) break;
		if (rsl_qsweep[tk_sindex] == 0) continue;
	  }
	  /* If data for this parm type exists in this toolkit sweep,
		   then move it into a rsl sweep. */
	  if (vs->tk.ncell[tk_sindex][pindex] > 0)
		{
			sindex++;
		  v->sweep[sindex] = SweepBuild(gvl1, vs, calibr, vindex, pindex,
																		tk_sindex);
			v->sweep[sindex]->h.sweep_num = sindex + 1;
			if (radar_verbose_flag)
			  fprintf(stderr, "  rsl_sweep[%02d]  elev=%4.1f  nrays=%d  cells/ray=%d\n", 
								v->sweep[sindex]->h.sweep_num-1, v->sweep[sindex]->h.elev,
								vs->tk.nray[tk_sindex], vs->tk.ncell[tk_sindex][pindex]);
		}
	} /* end for (tk_sindex=0;...*/
	return(v);
}

/*************************************************************/
/*                                                           */
/*                       parmIdentify                        */
/*                                                           */
/*************************************************************/
int parmIdentify(char *parmName)
/* Identify the parameter type stored in the L1B_1C_GV structure.
	 Upon success, return the corresponding RSL radar volume XX_INDEX
	 value.
	 Upon failure, return -1 .
*/
{
	int vindex;

	if (strcmp(parmName, "Z") == 0)
	{
		vindex = DZ_INDEX;   invf = DZ_INVF;   f = DZ_F;
	}
	else if (strcmp(parmName, "V") == 0)
	{
		vindex = VR_INDEX;   invf = VR_INVF;   f = VR_F;
	}
	else if (strcmp(parmName, "QCZ") == 0)
	{
		vindex = CZ_INDEX;   invf = CZ_INVF;   f = CZ_F;
	}
	else if (strcmp(parmName, "ZDR") == 0)
	{
		vindex = ZD_INDEX;   invf = ZD_INVF;   f = ZD_F;
	}
	else if (strcmp(parmName, "QCZDR") == 0)
	{
		vindex = CD_INDEX;   invf = CD_INVF;   f = CD_F;
	}
	else if (strcmp(parmName, "QCMZ") == 0)
	{
		vindex = MZ_INDEX;   invf = MZ_INVF;   f = MZ_F;
	}
	else if (strcmp(parmName, "QCMZDR") == 0)
	{
		vindex = MD_INDEX;   invf = MD_INVF;   f = MD_F;
	}
	else  /* Unknown */
	{
		return(-1);
	}

	return(vindex);
}

/*************************************************************/
/*                                                           */
/*                    Radar_headerFill                       */
/*                                                           */
/*************************************************************/
void Radar_headerFill(Radar *radar, L1B_1C_GV *gvl1)
{
	double x;
	
	radar->h.month = (int)gvl1->volDes.month;
	radar->h.day = (int)gvl1->volDes.day;
	radar->h.year = (int)gvl1->volDes.year;
	radar->h.hour = (int)gvl1->volDes.hour;
	radar->h.minute = (int)gvl1->volDes.minute;
	radar->h.sec = (float)gvl1->volDes.second;
	strncpy(radar->h.radar_type, "**", 48);       /*********/
	radar->h.nvolumes = MAX_RADAR_VOLUMES;
	radar->h.number = MISSING_VAL;
	strncpy(radar->h.name, gvl1->sensor.radarDesc.radarName, 7);
	strncpy(radar->h.radar_name, gvl1->sensor.sweepInfo[0].radarName, 7);

	/* Radar Latitude */
	x = fabs(gvl1->sensor.radarDesc.radarLat);
	radar->h.latd = (int)floor(x);
	x = (x - radar->h.latd) * 60.0;
	radar->h.latm = (int)floor(x);
	x = (x - radar->h.latm) * 60.0;
	radar->h.lats = (int)floor(x + 0.5);  /* round up */
	if (gvl1->sensor.radarDesc.radarLat < 0)
	{
	  radar->h.latd = -radar->h.latd;
		radar->h.latm = -radar->h.latm;
		radar->h.lats = -radar->h.lats;
	}
	
	/* Radar Longitude */
	x = fabs(gvl1->sensor.radarDesc.radarLon);
	radar->h.lond = (int)floor(x);
	x = (x - radar->h.lond) * 60.0;
	radar->h.lonm = (int)floor(x);
	x = (x - radar->h.lonm) * 60.0;
	radar->h.lons = (int)floor(x + 0.5);  /* round up */
	if (gvl1->sensor.radarDesc.radarLon < 0)
	{
	  radar->h.lond = -radar->h.lond;
		radar->h.lonm = -radar->h.lonm;
		radar->h.lons = -radar->h.lons;
	}

	radar->h.height = (int)(1000.0 * gvl1->sensor.radarDesc.radarAlt + 0.5);
	radar->h.spulse = MISSING_VAL;  /* ns */
	radar->h.lpulse = MISSING_VAL;  /* ns */
}

/*************************************************************/
/*                                                           */
/*                        RadarBuild                         */
/*                                                           */
/*************************************************************/
Radar *RadarBuild(L1B_1C_GV *gvl1, VosSize *vs, float zCal)
/* Creates and fills a RSL radar structure with data obtained
	 from the L1B_1C_GV structure.

	 If success, returns a pointer to the radar structure.
	 If failure, returns NULL.
*/
{
	Radar *radar;
	extern int rsl_qfield[];
	int pindex, vindex;
	
	if (radar_verbose_flag)
	{
	  fprintf(stderr, "\n****** Moving VOS from toolkit L1GV structure -> RSL structure ...\n");
	}
	/* Create a structure of type Radar */
	radar = (Radar *)RSL_new_radar(MAX_RADAR_VOLUMES);
	if (radar == NULL) 
	{
		perror("RadarBuild(): Error creating radar structure.\n");
		return(NULL);
	}
	/* Initialize the Radar_header values. */
	Radar_headerFill(radar, gvl1);

	/* Build each of the 'nparm' volumes of the radar structure. */
	for (pindex=0; pindex<vs->tk.nparm; pindex++)
	{
		/* Identify parameter type, so we know which RSL volume to load the data
			 into. */
		vindex = parmIdentify(gvl1->sensor.parm[pindex]->parmDesc.parmName);
		if (vindex < 0)
		{
			fprintf(stderr, 
				  "RadarBuild(): Unexpected parameter type: %s found in HDF file.\n",
							gvl1->sensor.parm[pindex]->parmDesc.parmName);
		}
		/* Don't build mask volumes. */
		else if ((vindex == MZ_INDEX) || (vindex == MD_INDEX)) continue;
		else if (rsl_qfield[vindex] == 0) /* Don't build unselected volumes. */
		{
		  if (radar_verbose_flag)
			{
			  fprintf(stderr, "Field %s not selected for retrieval from HDF file.\n",
								gvl1->sensor.parm[pindex]->parmDesc.parmName);
				if (vindex == CZ_INDEX)
			    fprintf(stderr, "Field 'DZ' unselected for retrieval from 1C-51 file.\n");
				else if (vindex == CD_INDEX)
			    fprintf(stderr, "Field 'ZD' unselected for retrieval from 1C-51 file.\n");
			}
		}  /* end else if (rsl_qfield[vindex] == 0) */
		else if (vindex == CZ_INDEX)  /* Handle CZ and DZ volumes. */
		{
			/* Build the RSL CZ volume. */
			radar->v[vindex] = VolumeBuild(gvl1, vs, zCal, vindex, pindex);
			/* If required, build a RSL DZ volume. */
			if (rsl_qfield[DZ_INDEX])
			{
				if (radar_verbose_flag)
			    fprintf(stderr, "Constructing reflectivity volume 'DZ'\n");
				radar->v[DZ_INDEX] = VolumeBuild(gvl1, vs, zCal, DZ_INDEX, pindex);
			}
		}  /* end if (vindex == CZ_INDEX) */
		else if (vindex == CD_INDEX)  /* Handle CD and ZD volumes. */
		{
			/* Build the RSL CD volume. */
			radar->v[vindex] = VolumeBuild(gvl1, vs, 0.0, vindex, pindex);
			/* If required, build a RSL ZD volume. */
			if (rsl_qfield[ZD_INDEX])
			{
				if (radar_verbose_flag)
			    fprintf(stderr, "Constructing reflectivity volume 'ZD'\n");
				radar->v[ZD_INDEX] = VolumeBuild(gvl1, vs, 0.0, ZD_INDEX, pindex);
			}
		}  /* end if (vindex == CD_INDEX) */
		else   /* Handle all 1B-51 fields. (DZ, ZD, VR) */
		{
		  radar->v[vindex] = VolumeBuild(gvl1, vs, 0.0, vindex, pindex);
		}
	}  /* end for (pindex=0; ...) */

	return(radar);
}

/*************************************************************/
/*                                                           */
/*                       commentsRead                        */
/*                                                           */
/*************************************************************/
int commentsRead(VosSize *vs, float *zCal, char *comments, int productID)
{
/* Parse the comments field of the 'L1B_1C_GV' structure.
	 Retrieve the number_of_cells/ray values for each parameter, and
	 store in the 'VosSize' structure.

	 Returns: OK if success.
	          <0 if failure.
*/
	char *spointer;
	char record[2][2048];  /* 2 records is maximum possible. */
	char parseString[1024];
	int nrecords, pindex, tk_sindex;
	float qcParm[NUMBER_QC_PARAMS];

	/* Construct a format string to read the records in the comments 
		 field. A logical record here is the block of ascii characters
		 which details the toolkit dimensions of one VOS.
		 For a 1B-51 file, there should be 1 such record.
		 For a 1C-51 file, there is one additional record for QC parms.*/

	strcpy(parseString, "");
	strcat(parseString, "%[^*] %*[*\n]");
	nrecords = 1;
	if (productID == TK_L1C_GV)
	{
	  strcat(parseString, "%[^\n]");
		nrecords++;
	}
	/* Read all records from the comments field into the record buffers. */
	if (sscanf(comments, parseString, record[0], record[1]) != nrecords)
	  goto quit;


	if (sscanf(record[0], "nSweep=%d", &vs->tk.nsweep) != 1) goto quit;

	strcpy(parseString, "nRay=%d\n");
	for (pindex=0; pindex<vs->tk.nparm; pindex++)
    strcat(parseString, "nCell_parm[%*d]=%d\n");

	spointer = record[0];
	for (tk_sindex=0; tk_sindex<vs->tk.nsweep; tk_sindex++)
	{
		spointer = strstr(spointer, "nRay=");
		if (sscanf(spointer, parseString, &vs->tk.nray[tk_sindex],
							 &vs->tk.ncell[tk_sindex][0], &vs->tk.ncell[tk_sindex][1],
							 &vs->tk.ncell[tk_sindex][2], &vs->tk.ncell[tk_sindex][3])
				!= vs->tk.nparm+1) goto quit;
		spointer = spointer + 5;
	}


	/* If 1C-51 file, read the QC parameters into the qcParm array. */
	if (productID == TK_L1C_GV)
	{
		if (sscanf(record[1], "-hThresh1 %f -hThresh2 %f -hThresh3 %f -zThresh0 %f -zThresh1 %f -zThresh2 %f -zThresh3 %f -hFreeze %f -dbzNoise %f -zCal %f", 
							 &qcParm[HTHRESH1], &qcParm[HTHRESH2], &qcParm[HTHRESH3],
							 &qcParm[ZTHRESH0], &qcParm[ZTHRESH1], &qcParm[ZTHRESH2], &qcParm[ZTHRESH3],
							 &qcParm[HFREEZE], &qcParm[DBZNOISE], &qcParm[ZCAL]) != NUMBER_QC_PARAMS)
		  goto quit;
		/* Print out the QC parameters we've just read in. */
		/*
		if (radar_verbose_flag)
		{
			fprintf(stderr, "\n****** Reading VOS QC Parameters from HDF file...\n");
			fprintf(stderr, "hThresh1: %.2f   hThresh2: %.2f   hThresh3: %.2f\n",
							qcParm[HTHRESH1], qcParm[HTHRESH2], qcParm[HTHRESH3]);
			fprintf(stderr, "zThresh0: %.2f  zThresh1: %.2f  zThresh2: %.2f  zThresh3: %.2f\n",
							qcParm[ZTHRESH0], qcParm[ZTHRESH1], qcParm[ZTHRESH2], qcParm[ZTHRESH3]);
			fprintf(stderr, "hFreeze: %.2f  dbzNoise: %.2f   zCal: %.2f\n\n", 
							qcParm[HFREEZE], qcParm[DBZNOISE], qcParm[ZCAL]);
		}
		*/
		if (qcParm[ZCAL] <= NOVAL_FLOAT) *zCal = 0.0;
		else *zCal = qcParm[ZCAL];
	} /* end if (productID == TK_L1C_GV) */

	return(OK);

 quit:
	if (radar_verbose_flag)
	  fprintf(stderr, "commentsRead(): Failure reading comments field\n");
	return(ABORT);
}

/*************************************************************/
/*                                                           */
/*                          GVL1Build                        */
/*                                                           */
/*************************************************************/
static L1B_1C_GV *GVL1Build(IO_HANDLE *granuleHandle, int vosNum,
														VosSize *vs)
{
/* Build a toolkit 'L1B_1C_GV' structure sized for the VOS we 
	 will later read in from the HDF file.
	 Returns:
	   gvl1 if success.
		 NULL if fails.
*/
	int ncell, pindex;
	L1B_1C_GV *gvl1;

	/* Using the toolkit, get the toolkit VOS dimensions from
		 the HDF file. Note that the toolkit dimensions are distinct
		 from the RSL VOS dimensions. */
	vs->tk.nparm = TKgetNparm(granuleHandle, vosNum);
/*  TK_FAIL is now defined as 1?????? */
	if (vs->tk.nparm <= 0)
	{
		fprintf(stderr, "GVL1Build(): TKgetNparm() failed.\n");
		return(NULL);
	}
	vs->tk.nsweep = TKgetNsweep(granuleHandle, vosNum);
	if (vs->tk.nsweep <= 0)
	{
		fprintf(stderr, "GVL1Build(): TKgetNsweep() failed.\n");
		return(NULL);
	}
	vs->rsl.maxNray = TKgetNray(granuleHandle, vosNum);
	if (vs->rsl.maxNray <= 0)
	{
		fprintf(stderr, "GVL1Build(): TKgetNray() failed.\n");
		return(NULL);
	}

	/* Allocate memory for a TSDIS 'L1B_1C_GV' structure. */
	gvl1 = (L1B_1C_GV *)TKnewGVL1();
	for (pindex=0; pindex<vs->tk.nparm; pindex++)
	{
		/* Allocate memory for a new parameter. */
		gvl1->sensor.parm[pindex] = (PARAMETER *)TKnewGVL1parm();
		ncell = TKgetNcell(granuleHandle, vosNum, pindex);
		/* Allocate memory for a 3D array to contain mask values and/or data. */
		if (granuleHandle->productID == TK_L1B_GV)
		{
		  gvl1->sensor.parm[pindex]->parmData2byte = 
			    (int16 ***)TKnewParmData2byte(vs->tk.nsweep, vs->rsl.maxNray, ncell);
		}
		else   /* 1C-51 */
		{
			/* Odd parameters contain data, even parameters contain masks. */
			if ((pindex % 2) == 0)  /* Mask? */
		    gvl1->sensor.parm[pindex]->parmData1byte = 
				   (int8 ***)TKnewParmData1byte(vs->tk.nsweep, vs->rsl.maxNray, ncell);
			else  /* data */
		    gvl1->sensor.parm[pindex]->parmData2byte = 
				   (int16 ***)TKnewParmData2byte(vs->tk.nsweep, vs->rsl.maxNray, ncell);
		}
	} /* end for (pindex=0; ... */
	
	return(gvl1);
}

/*************************************************************/
/*                                                           */
/*                   metaDataRead                            */
/*                                                           */
/*************************************************************/
int metaDataRead(Radar *radar, IO_HANDLE *granuleHandle)
{
	char buf[64];

	TKreadMetadataChar(granuleHandle, TK_RADAR_CITY, buf);
	strncpy(radar->h.city, buf, 14);
	TKreadMetadataChar(granuleHandle, TK_RADAR_STATE, buf);
	strncpy(radar->h.state, buf, 2);
	return(OK);
}

/*************************************************************/
/*                                                           */
/*                        hdfFileOpen                        */
/*                                                           */
/*************************************************************/
static int hdfFileOpen(char *infile, IO_HANDLE *granuleHandle, 
											 char *hdfFileName, int *vosNum)
{
/* Opens, if necessary, an HDF file. Checks that a VOS, not previously
	 retrieved, exists in the HDF file.

   Returns:
	   OK, if success.
		 <0, if failure.
*/
	char *product;
	int productID, nvos, status;

	/* If we presently have an open HDF file, check that it is the 
	   file requested; ie, 'infile'. If it's not, we first close
		 the open HDF file. */
	if (*vosNum != 0)
	{
		if (strcmp(hdfFileName, infile) != 0)
		{
			if (TKclose(granuleHandle) != TK_SUCCESS)
			{
			  fprintf(stderr, "hdfFileOpen(): *** TKclose() error\n");
				return(ABORT);
			}
			*vosNum = 0;
		}
	} /* end if (*vosNum != 0) */
	
	/* If first VOS of HDF file, we need first to open the file. */
	if (*vosNum == 0)
	{
		strncpy(hdfFileName, infile, TK_MAX_FILENAME-1);
		/* Get the desired product out of the HDF filename. */
		product = strrchr(hdfFileName, '/');
		if (product == NULL)
		  product = hdfFileName;
		else 
		  product = (char *)(product + 1);
		if (strncmp(product, "1B51", 4) == 0)
		  productID = TK_L1B_GV;
		else if (strncmp(product, "1C51", 4) == 0)
		  productID = TK_L1C_GV;
		else
		{
			fprintf(stderr, "hdfFileOpen(): Invalid HDF filename.\n");
			return(ABORT);
		}
		status = TKopen(hdfFileName, productID, TK_READ_ONLY, granuleHandle); 
		if (status != TK_SUCCESS)
		{
			fprintf(stderr, "hdfFileOpen(): *** TKopen() error\n");
			return(ABORT);
		}
	}  /* end if (*vosNum == 0) */
	
	/* Check if the requested VOS exists in the HDF granule. */
	nvos = (int)TKgetNvos(granuleHandle);
	if (nvos == 0)
	{
		if (radar_verbose_flag)
		  fprintf(stderr, "\nEmpty granule.\n");
		return(QUIT);
	}		
	else if (*vosNum+1 > nvos)
	{
		if (radar_verbose_flag)
		  fprintf(stderr, "\nAll VOSs read from HDF file: %s\n", hdfFileName);
		return(QUIT);
	}
	else if (nvos < 0)
	{
		fprintf(stderr, "hdfFileOpen():*** TKgetNvos() error\n");
		return(QUIT);
	}

	return(OK);
}

/*************************************************************/
/*                                                           */
/*                      RSL_hdf_to_radar                     */
/*                                                           */
/*************************************************************/
Radar *RSL_hdf_to_radar(char *infile)
{
	/* Reads one volume scan from a HDF file into a RSL radar structure.
		 It is envisioned that all VOSs will normally be retrieved,
		 one after the other, from an HDF file. Therefore, in the
		 interest of efficiency, this function is designed to keep
		 an HDF file open between VOS retrievals.

		 Returns:
		   - pointer to a filled radar structure, if success.
			 - NULL pointer, if failure.
	*/
	int pindex, status;
	float zCal=0.0;         /* Z Calibration constant. */
	VosSize vs;             /* VOS dimensions storage. */
	Radar *radar;           /* RSL structure for VOS storage. */
	L1B_1C_GV *gvl1;    /* TSDIS structure for VOS storage. */
	/* Following values must remain between invocations of this function. */
	static IO_HANDLE granuleHandle;
	static char hdfFileName[TK_MAX_FILENAME];
	static int vosNum=0;  /* No. of last VOS read from HDF file. (1...N) */
	
	/* Open, if necessary, an HDF file. */
	status = hdfFileOpen(infile, &granuleHandle, hdfFileName, &vosNum);
	if (status < 0) goto quit;

	/* Initialize the 'VosSize' structure. */
	memset(&vs, '\0', sizeof(VosSize));

	/* Build a toolkit 'L1B_1C_GV' structure correctly sized for the VOS we 
	 will next read in from the HDF file. */
	gvl1 = (L1B_1C_GV *)GVL1Build(&granuleHandle, vosNum, &vs);
	if (gvl1 == NULL) goto quit;
	
	/* Read VOS from HDF file into the toolkit L1B_1C_GV structure. */
	if (radar_verbose_flag)
	  fprintf(stderr, "\n\n***** Moving VOS from HDF file -> toolkit L1GV structure ...\n");
	status = TKreadL1GV(&granuleHandle, gvl1);
	if (status != TK_SUCCESS)
	{
		fprintf(stderr, "RSL_hdf_to_radar(): *** TKreadL1GV() error\n");
		goto quit;
	}
	
	if (radar_verbose_flag)
	{
		fprintf(stderr, "Input file: %s\n", hdfFileName);
		fprintf(stderr, "VOS date:   %.2d/%.2d/%d\n", gvl1->volDes.month, 
						gvl1->volDes.day, gvl1->volDes.year);
		fprintf(stderr, "VOS time:   %.2d:%.2d:%.2d\n", gvl1->volDes.hour, 
						gvl1->volDes.minute, gvl1->volDes.second);

		fprintf(stderr, "Granule VOS #: %d\n", gvl1->volDes.volNum);
		fprintf(stderr, "VOS Fields:\n");
		for (pindex=0; pindex<vs.tk.nparm ; pindex++)
		  fprintf(stderr, "  %d: %s\n", pindex+1,
							gvl1->sensor.parm[pindex]->parmDesc.parmDesc);
	}

	/* Scan thru the comments field of the 'gvl1' structure for
		 the number_of_cells/ray/sweep for each parameter, and store
		 the values in vs, so we can next build a correctly sized RSL
		 structure to contain the VOS.
	*/
	if (commentsRead(&vs, &zCal, gvl1->comments, granuleHandle.productID) < 0)
	  goto quit;

	/* Move VOS from L1B_1C_GV structure to radar structure. */
	radar = (Radar *)RadarBuild(gvl1, &vs, zCal);
	/* There are a couple metadata items needed for insertion into the
		 radar->header. */
	status = metaDataRead(radar, &granuleHandle);

	/* Free memory allocated to the toolkit GVL1 structure. */
	TKfreeGVL1(gvl1);
	if (radar == NULL) goto quit;

	vosNum++;
	return(radar);

 quit:
	if (status == QUIT)
	  TKclose(&granuleHandle);
	return(NULL);
}

#else
/*
 * Just declare and return something when we're told we don't have
 * TSDISTK nor HDF.  Do this because RSL_anyformat_to_radar references
 * this routine; linking won't fail because of no HDF.
 */

#include "rsl.h"
Radar *RSL_hdf_to_radar(char *infile)
{
  return NULL;
}
#endif
