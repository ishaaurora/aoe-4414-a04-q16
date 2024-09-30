# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
# _o represents the ECEF origin of SEZ frame and the others are the ecef position
# Converts ECEF vector components to SEZ
#  See "Fundamentals of Astrodynamics and Applications, Fourth Edition" by
#  David A. Vallado, pages 172-173
# Parameters:
# o_x_km: x ECEF origin of SEZ frame
# o_y_km y ECEF origin of SEZ frame
# o_z_km: z ECEF origin of SEZ frame
# x_km: ECEF position in x
# y_km: ECEF position in y
# z_km: ECEF position in z

# Output:
#  Prints the s e and z components in km
#
# Written by Isha Aurora
# Other contributors: None
#
# This work is licensed under CC BY-SA 4.0

# import Python modules
import math # math module
import sys  # argv

# "constants"
R_E_KM = 6378.1363
E_E    = 0.081819221456

# helper functions

## calculated denominator
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-(ecc**2)*(math.sin(lat_rad)**2))

# initialize script arguments
o_x_km = float('nan')  
o_y_km = float('nan')  
o_z_km = float('nan')  
x_km = float('nan')
y_km = float('nan')
z_km = float('nan')

# parse script arguments (always 1 more than the number of arguments)
if len(sys.argv) == 7:
    try:
        o_x_km = float(sys.argv[1])
        o_y_km = float(sys.argv[2])
        o_z_km = float(sys.argv[3])
        x_km = float(sys.argv[4])
        y_km = float(sys.argv[5])
        z_km = float(sys.argv[6])
    except ValueError:
        print("Error: o_x_km, o_y_km, o_z_km, x_km, y_km, and z_km must be numeric.")
        exit()
else:
    print('python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km')
    exit()

# write script below this line

#ecef to llh
# calculate longitude
lon_rad = math.atan2(o_y_km,o_x_km)
lon_deg = lon_rad*180.0/math.pi
# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
r_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')
# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
    denom = calc_denom(E_E,lat_rad)
    c_E = R_E_KM/denom
    prev_lat_rad = lat_rad
    lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
    count = count+1
# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E

#ecef_to_sez
#ecef vector from station to satellite
x_ecef_km = x_km - o_x_km
y_ecef_km = y_km - o_y_km
z_ecef_km = z_km - o_z_km

#get sez by doing rotation matrix
sin_lat_o = math.sin(lat_rad)
cos_lat_o = math.cos(lat_rad)
sin_lon_o = math.sin(lon_rad)
cos_lon_o = math.cos(lon_rad)

s_km = -z_ecef_km * cos_lat_o + x_ecef_km * cos_lon_o * sin_lat_o + y_ecef_km * sin_lat_o * sin_lon_o
e_km = y_ecef_km * cos_lon_o - x_ecef_km * sin_lon_o
z_km = x_ecef_km * cos_lat_o * cos_lon_o + z_ecef_km * sin_lat_o + y_ecef_km * cos_lat_o * sin_lon_o

#print SEZ coords
print(s_km)
print(e_km)
print(z_km)