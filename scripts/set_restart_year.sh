#!/bin/bash
#  Copyright 2020 Scott Wales
#
#  \author  Scott Wales <scott.wales@unimelb.edu.au>
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

# Sets the start year in each model component

source /etc/profile.d/modules.sh
module use /g/data/hh5/public/modules
module load conda/analysis3
module load nco

set -eu
trap "echo Error in set_restart_year.sh" ERR
export UMDIR=~access/umdir

# Load some helper scripts
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source $SCRIPTDIR/utils.sh

# Starting year
start_year=$1

# Set the restart year in the atmosphere and ice namelists
set_um_start_year $start_year
set_cice_start_year $start_year

# Most recent restart directory
payu_restart=$(ls -d ./archive/restart* | sort -t t -k 3 -n | tail -n 1)

echo "Setting restart year in ${payu_restart} to ${start_year}"

if [ ! -d $payu_restart/ocean ]; then
    echo "No restart directory"
    exit 1
fi

# Update ocean start time - read from a text file
cat > $payu_restart/ocean/ocean_solo.res << EOF
     3        (Calendar: no_calendar=0, thirty_day_months=1, julian=2, gregorian=3, noleap=4)
     1     1     1     0     0     0        Model start time:   year, month, day, hour, minute, second
   102     1     1     0     0     0        Current model time: year, month, day, hour, minute, second
EOF

# Update atmos start time - field in the restart file
python scripts/update_um_year.py $start_year $payu_restart/atmosphere/restart_dump.astart 2> /dev/null

cat > $payu_restart/atmosphere/um.res.yaml << EOF
end_date: $(printf %04d $start_year)-01-01 00:00:00
EOF

cat > $payu_restart/ice/restart_date.nml << EOF
&coupling
    inidate = 1010101
    init_date = 10101
/
EOF

# Clear ice step count
cat > $payu_restart/ice/cice_in.nml << EOF
&setup_nml
    days_per_year = 365
    year_init = 1
    istep0 = 0
    dt = 3600
    npt = 8760.0
    ndyn_dt = 1
    runtype = 'continue'
    ice_ic = 'default'
    restart = .true.
    restart_dir = './RESTART/'
    restart_file = 'iced'
    pointer_file = './RESTART/ice.restart_file'
    dumpfreq = 'y'
    dumpfreq_n = 1
    diagfreq = 24
    diag_type = 'file'
    diag_file = 'ice_diag.d'
    print_global = .false.
    print_points = .false.
    latpnt(1:2) = 90.0, -65.0
    lonpnt(1:2) = 0.0, -45.0
    dbug = .false.
    history_dir = './HISTORY/'
    write_ic = .false.
    incond_dir = './HISTORY/'
    incond_file = 'iceh_ic'
    history_file = 'iceh'
    history_format = 'nc'
    histfreq = 'd', 'm', 'x', 'x', 'x'
    histfreq_n = 1, 1, 1, 1, 1
    hist_avg = .true.
/

&grid_nml
    grid_format = 'nc'
    grid_type = 'tripole'
    grid_file = './INPUT/grid.nc'
    kmt_file = './INPUT/kmt.nc'
    kcatbound = 0
/

&domain_nml
    nprocs = 12
    processor_shape = 'slenderX1'
    distribution_type = 'cartesian'
    distribution_wght = 'latitude'
    ew_boundary_type = 'cyclic'
    ns_boundary_type = 'tripole'
/

&tracer_nml
    tr_iage = .false.
    restart_age = .false.
    tr_lvl = .false.
    restart_lvl = .false.
    tr_pond = .false.
    restart_pond = .false.
/

&ice_nml
    kitd = 1
    kdyn = 1
    ndte = 120
    kstrength = 1
    krdg_partic = 1
    krdg_redist = 1
    mu_rdg = 3
    advection = 'remap'
    heat_capacity = .false.
    conduct = 'bubbly'
    shortwave = 'default'
    albedo_type = 'default'
    albicev = 0.86
    albicei = 0.44
    albsnowv = 0.98
    albsnowi = 0.7
    ahmax = 0.1
    snowpatch = 0.01
    dt_mlt = 1.0
    dalb_mlt = -0.02
    awtvdr = 0.00318
    awtidr = 0.00182
    awtvdf = 0.63282
    awtidf = 0.36218
    r_ice = 0.0
    r_pnd = 0.0
    r_snw = 0.0
    atmbndy = 'default'
    fyear_init = 1
    ycycle = 1
    atm_data_format = 'nc'
    atm_data_type = 'default'
    atm_data_dir = 'unknown_atm_data_dir'
    calc_strair = .false.
    calc_tsfc = .false.
    precip_units = 'mks'
    tfrzpt = 'linear_S'
    tocnfrz = -1.8
    ustar_min = 0.0005
    cosw = 0.96
    sinw = 0.28
    dragio = 0.00536
    chio = 0.004
    iceruf = 0.0005
    update_ocn_f = .true.
    oceanmixed_ice = .false.
    ocn_data_format = 'nc'
    sss_data_type = 'default'
    sst_data_type = 'default'
    ocn_data_dir = 'unknown_ocn_data_dir'
    oceanmixed_file = 'unknown_oceanmixed_file'
    restore_sst = .false.
    trestore = 0
    restore_ice = .false.
/

&icefields_nml
    f_tmask = .true.
    f_tarea = .true.
    f_uarea = .true.
    f_dxt = .false.
    f_dyt = .false.
    f_dxu = .false.
    f_dyu = .false.
    f_htn = .false.
    f_hte = .false.
    f_angle = .true.
    f_anglet = .true.
    f_ncat = .true.
    f_vgrdi = .true.
    f_vgrds = .true.
    f_bounds = .false.
    f_hi = 'm'
    f_hs = 'm'
    f_tsfc = 'x'
    f_aice = 'm'
    f_uvel = 'm'
    f_vvel = 'm'
    f_fswdn = 'm'
    f_flwdn = 'm'
    f_snow = 'x'
    f_snow_ai = 'm'
    f_rain = 'x'
    f_rain_ai = 'm'
    f_sst = 'm'
    f_sss = 'm'
    f_uocn = 'x'
    f_vocn = 'x'
    f_frzmlt = 'x'
    f_fswfac = 'x'
    f_fswabs = 'x'
    f_fswabs_ai = 'm'
    f_albsni = 'x'
    f_alvdr = 'x'
    f_alidr = 'x'
    f_albice = 'x'
    f_albsno = 'x'
    f_albpnd = 'x'
    f_coszen = 'x'
    f_flat = 'x'
    f_flat_ai = 'm'
    f_fsens = 'x'
    f_fsens_ai = 'm'
    f_flwup = 'x'
    f_flwup_ai = 'm'
    f_evap = 'x'
    f_evap_ai = 'm'
    f_tair = 'x'
    f_tref = 'x'
    f_qref = 'x'
    f_congel = 'x'
    f_frazil = 'x'
    f_snoice = 'x'
    f_meltt = 'x'
    f_meltb = 'x'
    f_meltl = 'x'
    f_fresh = 'x'
    f_fresh_ai = 'm'
    f_fsalt = 'x'
    f_fsalt_ai = 'x'
    f_fhocn = 'x'
    f_fhocn_ai = 'm'
    f_fswthru = 'x'
    f_fswthru_ai = 'x'
    f_fsurf_ai = 'm'
    f_fcondtop_ai = 'm'
    f_fmeltt_ai = 'm'
    f_strairx = 'm'
    f_strairy = 'm'
    f_strtltx = 'm'
    f_strtlty = 'm'
    f_strcorx = 'm'
    f_strcory = 'm'
    f_strocnx = 'm'
    f_strocny = 'm'
    f_strintx = 'm'
    f_strinty = 'm'
    f_strength = 'm'
    f_divu = 'm'
    f_shear = 'm'
    f_sig1 = 'x'
    f_sig2 = 'x'
    f_dvidtt = 'x'
    f_dvidtd = 'x'
    f_daidtt = 'x'
    f_daidtd = 'x'
    f_mlt_onset = 'm'
    f_frz_onset = 'm'
    f_dardg1dt = 'x'
    f_dardg2dt = 'm'
    f_dvirdgdt = 'x'
    f_opening = 'm'
    f_hisnap = 'x'
    f_aisnap = 'x'
    f_trsig = 'x'
    f_icepresent = 'm'
    f_iage = 'x'
    f_alvl = 'x'
    f_vlvl = 'x'
    f_ardg = 'x'
    f_vrdg = 'x'
    f_aicen = 'm'
    f_vicen = 'm'
    f_tinz = 'x'
    f_tsnz = 'x'
    f_fsurfn_ai = 'm'
    f_fcondtopn_ai = 'm'
    f_fmelttn_ai = 'm'
    f_flatn_ai = 'm'
    f_apondn = 'x'
/
EOF

# Get the number of seconds since the run start date in ice/input_ice.nml
# init_date is the initial date of the experiment
# inidate is the date of the current run
runtime0=$(python <<EOF
import f90nml
import sys
import cftime
t0 = f90nml.read('ice/input_ice.nml')['coupling']['init_date']
y = t0//10000; m = (t0//100)%100; d = t0%100 
t0 = cftime.datetime(y,m,d,calendar='proleptic_gregorian')
t1 = cftime.datetime($start_year,1,1,calendar='proleptic_gregorian')
diff = t1 - t0
print(f"ice init_date {t0}, runtime0 {diff}", file=sys.stderr)
print(int(diff.total_seconds()))
EOF
)

# Put this in the coupling namelist - used by Payu to generate Oasis namcouple file
cat > $payu_restart/ice/input_ice.nml << EOF
&coupling
runtime0=$runtime0
runtime=0
/
EOF

# Set the date in the cice netcdf file
ncatted -a units,time,o,c,"seconds since ${start_year}-01-01 00:00:00" $payu_restart/ice/mice.nc

# Seconds between init_date and inidate
secs_realyr=0

ice_restart=$(ls $payu_restart/ice/iced.*0101)
mv $ice_restart ${ice_restart}.orig

# Set the date in the cice binary restart file
scripts/cicedumpdatemodify.py -i ${ice_restart}.orig -o $payu_restart/ice/iced.${start_year}0101 --istep0=0 --time=${secs_realyr}. --time_forc=0.

rm ${ice_restart}.orig
