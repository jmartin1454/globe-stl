#!/usr/bin/env python3

exaggeration   = 20
header         = ('Dry Earth %s-times-exaggerated elevation model by CMG Lee using NASA data.'%(exaggeration))
path_png_alt   = 'Earth_dry_elevation.png' ## 1-channel equirectangular PNG
luma_datum     = 141                       ## image intensity level (of 0-255) of datum
radius_datum   = 6378.137                  ## mean radius of zero level in km
f_wgs84        = 1 / 298.257223563         ## WGS84 flattening factor
km_per_luma    = (10.994 + 8.848) / 255 * exaggeration ## min and max elevations in km
scale          = 1e-2                      ## overall scale of model in km^-1
lat_offset     = 5.0 / 8                   ## rotation around planet axis in revolutions
n_division     = 800                       ## each cubic face divided into n_division^2 squares
# original value of n_division is 200 (Jeff)

class Png:
 def __init__(self, path):
  (self.width, self.height, self.pixels, self.metadatas) = png.Reader(path).read_flat()
 def __str__(self): return str((self.width, self.height, len(self.pixels), self.metadatas))

import time, re, math, struct, png
time.start = time.time()
def log(string): print('%6.3fs\t%s' % (time.time() - time.start, string))
def fmt(string): ## string.format(**vars()) using tags {expression!format} by CMG Lee
 def f(tag): i_sep = tag.rfind('!'); return (re.sub('\.0+$', '', str(eval(tag[1:-1])))
  if (i_sep < 0) else ('{:%s}' % tag[i_sep + 1:-1]).format(eval(tag[1:i_sep])))
 return (re.sub(r'(?<!{){[^{}]+}', lambda m:f(m.group()), string)
         .replace('{{', '{').replace('}}', '}'))
def append(obj, string): return obj.append(fmt(string))
def tabbify(cellss, separator='|'):
 cellpadss = [list(rows) + [''] * (len(max(cellss, key=len)) - len(rows)) for rows in cellss]
 fmts = ['%%%ds' % (max([len(str(cell)) for cell in cols])) for cols in zip(*cellpadss)]
 return '\n'.join([separator.join(fmts) % tuple(rows) for rows in cellpadss])
def hex_rgb(colour): ## convert [#]RGB to #RRGGBB and [#]RRGGBB to #RRGGBB
 return '#%s' % (colour if len(colour) > 4 else ''.join([c * 2 for c in colour])).lstrip('#')
def viscam_colour(colour):
 colour_hex      = hex_rgb(colour)
 colour_top5bits = [int(colour_hex[i:i+2], 16) >> 3 for i in range(1,7,2)]
 return (1 << 15) + (colour_top5bits[0] << 10) + (colour_top5bits[1] << 5) + colour_top5bits[2]
def roundm(x, multiple=1):
 if   (isinstance(x, tuple)): return tuple(roundm(list(x), multiple))
 elif (isinstance(x, list )): return [roundm(x_i, multiple) for x_i in x]
 else: return int(math.floor(float(x) / multiple + 0.5)) * multiple
def average(xs): return None if (len(xs) == 0) else float(sum(xs)) / len(xs)
def flatten(lss): return [l for ls in lss for l in ls]
def rotate(facetss, degs): ## around x then y then z axes
 (deg_x,deg_y,deg_z) = degs
 (sin_x,cos_x) = (math.sin(math.radians(deg_x)), math.cos(math.radians(deg_x)))
 (sin_y,cos_y) = (math.sin(math.radians(deg_y)), math.cos(math.radians(deg_y)))
 (sin_z,cos_z) = (math.sin(math.radians(deg_z)), math.cos(math.radians(deg_z)))
 facet_rotatess = []
 for facets in facetss:
  facet_rotates = []
  for i_point in range(4):
   (x,y,z) = [facets[3 * i_point + i_xyz] for i_xyz in range(3)]
   if (x is None or y is None or z is None): facet_rotates += [x,y,z]

   else:
    (y,z) = (y * cos_x - z * sin_x, y * sin_x + z * cos_x) ## rotate about x
    (x,z) = (x * cos_y + z * sin_y,-x * sin_y + z * cos_y) ## rotate about y
    (x,y) = (x * cos_z - y * sin_z, x * sin_z + y * cos_z) ## rotate about z
    facet_rotates += [round(value, 9) for value in [x,y,z]]
  facet_rotatess.append(facet_rotates)
 return facet_rotatess
def translate(facetss, ds): ## ds = (dx,dy,dz)
 return [facets[:3] + [facets[3 * i_point + i_xyz] + ds[i_xyz]
                       for i_point in range(1,4) for i_xyz in range(3)]  for facets in facetss]
def flip(facetss): return [facets[:3]+facets[6:9]+facets[3:6]+facets[9:] for facets in facetss]

def cube_xyz_to_sphere_xyz(cube_xyzs):
 (x,y,z)                         = [float(xyz) for xyz in cube_xyzs]
 (x_squared,y_squared,z_squared) = (x * x,y * y,z * z)
 return (x * (1 - (y_squared + z_squared) / 2 + y_squared * z_squared / 3) ** 0.5,
         y * (1 - (x_squared + z_squared) / 2 + x_squared * z_squared / 3) ** 0.5,
         z * (1 - (y_squared + x_squared) / 2 + y_squared * x_squared / 3) ** 0.5)
def xyz_to_lla(xyzs):
 (x,y,z) = xyzs
 alt     = (x * x + y * y + z * z) ** 0.5
 lon     = math.atan2(y, x)
 lat     = math.asin(z / alt)
 return (lat,lon,alt)
deg_90 = math.pi / 2
def find_alt(lat_lons, altss):
  (lat,lon) = lat_lons
  #print(lat,lon) # Jeff's print statement
  if   (lat ==  deg_90): alt = average(altss[ 0])
  elif (lat == -deg_90): alt = average(altss[-1])
  else:
   (width,height) = (len(altss[0]),len(altss))
   x              = (0.5 + lon / (deg_90 * 4) + lat_offset) * width
   y              = (0.5 - lat / (deg_90 * 2)             ) * height
   (x_int,y_int)  = (int(x)   , int(y)   )
   (x_dec,y_dec)  = (x - x_int, y - y_int)
   (x0,x1)        = (x_int % width , (x_int + 1) % width )
   (y0,y1)        = (y_int % height, (y_int + 1) % height)
   #print(width,height,x,y,x_int,y_int,x0,x1,y0,y1) # Jeff's print statement
   # Jeff's modified version:
   #alt            = ((altss[y0][x0] * (1 - y_dec) + altss[y1][x0] * y_dec) * (1 - x_dec) +
   #                  (altss[y0][x1] * (1 - y_dec) + altss[y1][x1] * y_dec) *      x_dec)
   #print(alt,x_dec,y_dec)
   # original version:
   alt            = ((altss[y0][x0] * (1 - x_dec) + altss[y1][x0] * x_dec) * (1 - y_dec) +
                     (altss[y0][x1] * (1 - x_dec) + altss[y1][x1] * x_dec) *      y_dec)
   # above line looks like a mistake, to me (Jeff)
   # print(map(math.degrees, lat_lons), y,x, alt)
   #print(alt,x_dec,y_dec,altss[y0][x0],altss[y1][x0],altss[y0][x1],altss[y1][x1])
  return alt
def radius_wgs84(lat):
 if (lat in radius_wgs84.cachess): return radius_wgs84.cachess[lat]
 (sin_lat, cos_lat)        = (math.sin(lat), math.cos(lat))
 ff                        = (1 - f_wgs84) ** 2
 c                         = 1 / (cos_lat ** 2 + ff * sin_lat ** 2) ** 0.5
 s                         = c * ff
 radius_c_s_s              = (radius_datum * c, radius_datum * s)
 radius_wgs84.cachess[lat] = radius_c_s_s
 return radius_c_s_s
radius_wgs84.cachess = {}
def lla_to_sphere_xyz(llas):
 (lat,lon,alt)        = llas
 (sin_lat,sin_lon)    = (math.sin(lat),math.sin(lon))
 (cos_lat,cos_lon)    = (math.cos(lat),math.cos(lon))
 (radius_c, radius_s) = [(c_s_radius + alt * km_per_luma) * scale
                         for c_s_radius in radius_wgs84(lat)]
 return (radius_c * cos_lat * cos_lon,radius_c * cos_lat * sin_lon,radius_s * sin_lat)
def xyz_alt_to_xyza(xyzs, altss):
 (lat,lon,alt) = xyz_to_lla(xyzs)
 alt           = find_alt((lat,lon), altss)
 lla_alts      = [list(lla_to_sphere_xyz((lat,lon,alt))), alt]
 return lla_alts

log("Read elevation data")
png_alt = Png(path_png_alt)
if (png_alt.metadatas['planes'] != 1): print("%s not 1-channel PNG" % (path_png_alt)); sys.exit(1)
log(png_alt)
altss = [[png_alt.pixels[png_alt.width * y + x] - luma_datum
          for x in range(png_alt.width)] for y in range(png_alt.height)] ## altss[y][x]

log("Find vertices")
k       = 2.0 / n_division
range_k = range(n_division + 1)
face_vertex_llassss = [ ## [0=top][i_y][i_x][xyz,alt]
 [[xyz_alt_to_xyza((x*k-1,y*k-1,    1), altss) for y in range_k] for x in range_k],
 [[xyz_alt_to_xyza((x*k-1,   -1,y*k-1), altss) for y in range_k] for x in range_k],
 [[xyz_alt_to_xyza((    1,x*k-1,y*k-1), altss) for y in range_k] for x in range_k],
 [[xyz_alt_to_xyza((y*k-1,x*k-1,   -1), altss) for y in range_k] for x in range_k],
 [[xyz_alt_to_xyza((y*k-1,    1,x*k-1), altss) for y in range_k] for x in range_k],
 [[xyz_alt_to_xyza((   -1,y*k-1,x*k-1), altss) for y in range_k] for x in range_k],
]

log("Add facets") ## cube xyz -> ll(a) -> image xy -> a -> sphere xyz
facetss = []
for (i_face,face_vertex_llasss) in enumerate(face_vertex_llassss):
 for  v in range(n_division):
  for u in range(n_division):
   (xyz00, alt00) = face_vertex_llasss[v    ][u    ]
   (xyz01, alt01) = face_vertex_llasss[v    ][u + 1]
   (xyz10, alt10) = face_vertex_llasss[v + 1][u    ]
   (xyz11, alt11) = face_vertex_llasss[v + 1][u + 1]
   (xyz_m, alt_m) = xyz_alt_to_xyza([average(xyzs) for xyzs in zip(*(xyz00,xyz01,xyz10,xyz11))],
                                    altss)
   if (alt_m > max(alt00,alt01,alt10,alt11) or alt_m < min(alt00,alt01,alt10,alt11)):
    facetss.append([None,0,0] + xyz_m + xyz00 + xyz10)
    facetss.append([None,0,0] + xyz_m + xyz10 + xyz11)
    facetss.append([None,0,0] + xyz_m + xyz11 + xyz01)
    facetss.append([None,0,0] + xyz_m + xyz01 + xyz00)
   else:
    if (abs(alt00 - alt11) < abs(alt01 - alt10)):
     facetss.append([None,0,0] + xyz00 + xyz10 + xyz11)
     facetss.append([None,0,0] + xyz11 + xyz01 + xyz00)
    else:
     facetss.append([None,0,0] + xyz10 + xyz11 + xyz01)
     facetss.append([None,0,0] + xyz01 + xyz00 + xyz10)

log("Calculate normals")
for facets in facetss:
 if (facets[0] is None or facets[1] is None or facets[2] is None):
  us      = [facets[i_xyz + 9] - facets[i_xyz + 6] for i_xyz in range(3)]
  vs      = [facets[i_xyz + 6] - facets[i_xyz + 3] for i_xyz in range(3)]
  normals = [us[1]*vs[2] - us[2]*vs[1], us[2]*vs[0] - us[0]*vs[2], us[0]*vs[1] - us[1]*vs[0]]
  normal_length = sum([component * component for component in normals]) ** 0.5
  facets[:3] = [-round(component / normal_length, 10) for component in normals]

# log(tabbify([['N%s'  % (xyz   )                   for xyz in list('xyz')] +
#              ['%s%d' % (xyz, n) for n in range(3) for xyz in list('XYZ')] + ['RGB']] + facetss))

log("Compile STL")

# Jeff rewrote the following so that very large file sizes are
# possible without storing everything in memory

f_out=open(__file__[:__file__.rfind('.')] + '-parts.stl', 'wb')

the_start=[[('STL\n\n%-73s\n\n' % (header[:73])).encode('utf-8'), struct.pack('<L',len(facetss))]]

outss = the_start
out   = b''.join([bytes(out) for outs in outss for out in outs])
f_out.write(out)

log("the start")

for facets in facetss:
 part1=[struct.pack('<f',float(value)) for value in facets[:12]]
 part2=[struct.pack('<H',0 if (len(facets) <= 12) else
                    viscam_colour(facets[12]))]
 outss=[part1+part2]
 out   = b''.join([bytes(out) for outs in outss for out in outs])
 f_out.write(out)

log("#bytes:%d\t#facets:%d\ttitle:\"%-73s\"" % (len(out), len(facetss), header[:73]))
