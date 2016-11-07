function LatLon(lat, lon, height) {
  if (arguments.length < 3) height = 0;
  this.lat = lat;
  this.lon = lon;
  this.height = height;
}
/*
 * convert geodesic co-ordinates to OS grid reference
 */
function LatLongToOSGrid(p) {
  var lat = p.lat.toRad(), lon = p.lon.toRad();

  var a = 6377563.396, b = 6356256.910;          // Airy 1830 major & minor semi-axes
  var F0 = 0.9996012717;                         // NatGrid scale factor on central meridian
  var lat0 = (49).toRad(), lon0 = (-2).toRad();  // NatGrid true origin
  var N0 = -100000, E0 = 400000;                 // northing & easting of true origin, metres
  var e2 = 1 - (b*b)/(a*a);                      // eccentricity squared
  var n = (a-b)/(a+b), n2 = n*n, n3 = n*n*n;

  var cosLat = Math.cos(lat), sinLat = Math.sin(lat);
  var nu = a*F0/Math.sqrt(1-e2*sinLat*sinLat);              // transverse radius of curvature
  var rho = a*F0*(1-e2)/Math.pow(1-e2*sinLat*sinLat, 1.5);  // meridional radius of curvature
  var eta2 = nu/rho-1;

  var Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (lat-lat0);
  var Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(lat-lat0) * Math.cos(lat+lat0);
  var Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(lat-lat0)) * Math.cos(2*(lat+lat0));
  var Md = (35/24)*n3 * Math.sin(3*(lat-lat0)) * Math.cos(3*(lat+lat0));
  var M = b * F0 * (Ma - Mb + Mc - Md);              // meridional arc

  var cos3lat = cosLat*cosLat*cosLat;
  var cos5lat = cos3lat*cosLat*cosLat;
  var tan2lat = Math.tan(lat)*Math.tan(lat);
  var tan4lat = tan2lat*tan2lat;

  var I = M + N0;
  var II = (nu/2)*sinLat*cosLat;
  var III = (nu/24)*sinLat*cos3lat*(5-tan2lat+9*eta2);
  var IIIA = (nu/720)*sinLat*cos5lat*(61-58*tan2lat+tan4lat);
  var IV = nu*cosLat;
  var V = (nu/6)*cos3lat*(nu/rho-tan2lat);
  var VI = (nu/120) * cos5lat * (5 - 18*tan2lat + tan4lat + 14*eta2 - 58*tan2lat*eta2);

  var dLon = lon-lon0;
  var dLon2 = dLon*dLon, dLon3 = dLon2*dLon, dLon4 = dLon3*dLon, dLon5 = dLon4*dLon, dLon6 = dLon5*dLon;

  var N = I + II*dLon2 + III*dLon4 + IIIA*dLon6;
  var E = E0 + IV*dLon + V*dLon3 + VI*dLon5;

  return gridrefNumToLet(E, N, 6);
}


/*
 * convert OS grid reference to geodesic co-ordinates
 */
/*
function OSGridToLatLong(gridRef) {
  var gr = gridrefLetToNum(gridRef);
  var E = gr[0], N = gr[1];

  var a = 6377563.396, b = 6356256.910;              // Airy 1830 major & minor semi-axes
  var F0 = 0.9996012717;                             // NatGrid scale factor on central meridian
  var lat0 = 49*Math.PI/180, lon0 = -2*Math.PI/180;  // NatGrid true origin
  var N0 = -100000, E0 = 400000;                     // northing & easting of true origin, metres
  var e2 = 1 - (b*b)/(a*a);                          // eccentricity squared
  var n = (a-b)/(a+b), n2 = n*n, n3 = n*n*n;

  var lat=lat0, M=0;
  do {
    lat = (N-N0-M)/(a*F0) + lat;

    var Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (lat-lat0);
    var Mb = (3*n + 3*n*n + (21/8)*n3) * Math.sin(lat-lat0) * Math.cos(lat+lat0);
    var Mc = ((15/8)*n2 + (15/8)*n3) * Math.sin(2*(lat-lat0)) * Math.cos(2*(lat+lat0));
    var Md = (35/24)*n3 * Math.sin(3*(lat-lat0)) * Math.cos(3*(lat+lat0));
    M = b * F0 * (Ma - Mb + Mc - Md);                // meridional arc

  } while (N-N0-M >= 0.00001);  // ie until < 0.01mm

  var cosLat = Math.cos(lat), sinLat = Math.sin(lat);
  var nu = a*F0/Math.sqrt(1-e2*sinLat*sinLat);              // transverse radius of curvature
  var rho = a*F0*(1-e2)/Math.pow(1-e2*sinLat*sinLat, 1.5);  // meridional radius of curvature
  var eta2 = nu/rho-1;

  var tanLat = Math.tan(lat);
  var tan2lat = tanLat*tanLat, tan4lat = tan2lat*tan2lat, tan6lat = tan4lat*tan2lat;
  var secLat = 1/cosLat;
  var nu3 = nu*nu*nu, nu5 = nu3*nu*nu, nu7 = nu5*nu*nu;
  var VII = tanLat/(2*rho*nu);
  var VIII = tanLat/(24*rho*nu3)*(5+3*tan2lat+eta2-9*tan2lat*eta2);
  var IX = tanLat/(720*rho*nu5)*(61+90*tan2lat+45*tan4lat);
  var X = secLat/nu;
  var XI = secLat/(6*nu3)*(nu/rho+2*tan2lat);
  var XII = secLat/(120*nu5)*(5+28*tan2lat+24*tan4lat);
  var XIIA = secLat/(5040*nu7)*(61+662*tan2lat+1320*tan4lat+720*tan6lat);

  var dE = (E-E0), dE2 = dE*dE, dE3 = dE2*dE, dE4 = dE2*dE2, dE5 = dE3*dE2, dE6 = dE4*dE2, dE7 = dE5*dE2;
  lat = lat - VII*dE2 + VIII*dE4 - IX*dE6;
  var lon = lon0 + X*dE - XI*dE3 + XII*dE5 - XIIA*dE7;

  return new LatLon(lat.toDeg(), lon.toDeg());
}
*/

/*
 * convert standard grid reference ('SU387148') to fully numeric ref ([438700,114800])
 *   returned co-ordinates are in metres, centred on grid square for conversion to lat/long
 *
 *   note that northern-most grid squares will give 7-digit northings
 *   no error-checking is done on gridref (bad input will give bad results or NaN)
 */
/*
function gridrefLetToNum(gridref) {
  // get numeric values of letter references, mapping A->0, B->1, C->2, etc:
  var l1 = gridref.toUpperCase().charCodeAt(0) - 'A'.charCodeAt(0);
  var l2 = gridref.toUpperCase().charCodeAt(1) - 'A'.charCodeAt(0);
  // shuffle down letters after 'I' since 'I' is not used in grid:
  if (l1 > 7) l1--;
  if (l2 > 7) l2--;

  // convert grid letters into 100km-square indexes from false origin (grid square SV):
  var e = ((l1-2)%5)*5 + (l2%5);
  var n = (19-Math.floor(l1/5)*5) - Math.floor(l2/5);

  // skip grid letters to get numeric part of ref, stripping any spaces:
  gridref = gridref.slice(2).replace(/ /g,'');

  // append numeric part of references to grid index:
  e += gridref.slice(0, gridref.length/2);
  n += gridref.slice(gridref.length/2);

  // normalise to 1m grid, rounding up to centre of grid square:
  switch (gridref.length) {
    case 6:e += '50';n += '50';break;
    case 8:e += '5';n += '5';break;
    // 10-digit refs are already 1m
  }
//    var gridStringE = e.toString();
//    var gridStringN = n.toString();
//    return gridStringE + gridStringN;
  return [e, n];

}
*/

/*
 * convert numeric grid reference (in metres) to standard-form grid ref
 */
function gridrefNumToLet(e, n, digits) {
  // get the 100km-grid indices
  var e100k = Math.floor(e/100000), n100k = Math.floor(n/100000);

  if (e100k<0 || e100k>6 || n100k<0 || n100k>12) return '';

  // translate those into numeric equivalents of the grid letters
  var l1 = (19-n100k) - (19-n100k)%5 + Math.floor((e100k+10)/5);
  var l2 = (19-n100k)*5%25 + e100k%5;

  // compensate for skipped 'I' and build grid letter-pairs
  if (l1 > 7) l1++;
  if (l2 > 7) l2++;
  var letPair = String.fromCharCode(l1+'A'.charCodeAt(0), l2+'A'.charCodeAt(0));

  // strip 100km-grid indices from easting & northing, and reduce precision
  e = Math.floor((e%100000)/Math.pow(10,5-digits/2));
  n = Math.floor((n%100000)/Math.pow(10,5-digits/2));

  var gridRef = letPair + e.padLZ(digits/2) + n.padLZ(digits/2);

  return gridRef;
}

/*
 * pad a number with sufficient leading zeros to make it w chars wide
 */
Number.prototype.padLZ = function(w) {
  var n = this.toString();
  var l = n.length;
  for (var i=0; i<(w-l); i++) {
      n = '0' + n;
  }
  return n;
}

String.prototype.parseDeg = function() {
  if (!isNaN(this)) return Number(this);                 // signed decimal degrees without NSEW

  var degLL = this.replace(/^-/,'').replace(/[NSEW]/i,'');  // strip off any sign or compass dir'n
  var dms = degLL.split(/[^0-9.]+/);                     // split out separate d/m/s
  for (var i in dms) if (dms[i]=='') dms.splice(i,1);    // remove empty elements (see note below)
  switch (dms.length) {                                  // convert to decimal degrees...
    case 3:                                              // interpret 3-part result as d/m/s
      var deg = dms[0]/1 + dms[1]/60 + dms[2]/3600;break;
    case 2:                                              // interpret 2-part result as d/m
      var deg = dms[0]/1 + dms[1]/60;break;
    case 1:                                              // decimal or non-separated dddmmss
      if (/[NS]/i.test(this)) degLL = '0' + degLL;       // - normalise N/S to 3-digit degrees
      var deg = dms[0].slice(0,3)/1 + dms[0].slice(3,5)/60 + dms[0].slice(5)/3600;break;
    default:return NaN;
  }
  if (/^-/.test(this) || /[WS]/i.test(this)) deg = -deg; // take '-', west and south as -ve
  return deg;
}
// note: whitespace at start/end will split() into empty elements (except in IE)


// extend Number object with methods for converting degrees/radians

Number.prototype.toRad = function() {  // convert degrees to radians
  return this * Math.PI / 180;
}

Number.prototype.toDeg = function() {  // convert radians to degrees (signed)
  return this * 180 / Math.PI;
}

// ---- the following are duplicated from latlong-convert-coords.html ---- //

// ellipse parameters
var e = {WGS84:    {a: 6378137,     b: 6356752.3142, f: 1/298.257223563},
          Airy1830: {a: 6377563.396, b: 6356256.910,  f: 1/299.3249646}};

// helmert transform parameters
var h = {WGS84toOSGB36: {tx: -446.448,  ty:  125.157,   tz: -542.060,   // m
                           rx:   -0.1502, ry:   -0.2470,  rz:   -0.8421,  // sec
                           s:    20.4894},                               // ppm
          OSGB36toWGS84: {tx:  446.448,  ty: -125.157,   tz:  542.060,
                           rx:    0.1502, ry:    0.2470,  rz:    0.8421,
                           s:   -20.4894}};

function convertWGS84toOSGB36(p1) {
  var p2 = convert(p1, e.WGS84, h.WGS84toOSGB36, e.Airy1830);
  return p2;
}

function convert(p, e1, t, e2) {
  // -- convert polar to cartesian coordinates (using ellipse 1)

  p1 = new LatLon(p.lat, p.lon, p.height);  // to avoid modifying passed param
  p1.lat = p.lat.toRad();p1.lon = p.lon.toRad();

  var a = e1.a, b = e1.b;

  var sinPhi = Math.sin(p1.lat), cosPhi = Math.cos(p1.lat);
  var sinLambda = Math.sin(p1.lon), cosLambda = Math.cos(p1.lon);
  var H = p1.height;

  var eSq = (a*a - b*b) / (a*a);
  var nu = a / Math.sqrt(1 - eSq*sinPhi*sinPhi);

  var x1 = (nu+H) * cosPhi * cosLambda;
  var y1 = (nu+H) * cosPhi * sinLambda;
  var z1 = ((1-eSq)*nu + H) * sinPhi;


  // -- apply helmert transform using appropriate params

  var tx = t.tx, ty = t.ty, tz = t.tz;
  var rx = t.rx/3600 * Math.PI/180;  // normalise seconds to radians
  var ry = t.ry/3600 * Math.PI/180;
  var rz = t.rz/3600 * Math.PI/180;
  var s1 = t.s/1e6 + 1;              // normalise ppm to (s+1)

  // apply transform
  var x2 = tx + x1*s1 - y1*rz + z1*ry;
  var y2 = ty + x1*rz + y1*s1 - z1*rx;
  var z2 = tz - x1*ry + y1*rx + z1*s1;


  // -- convert cartesian to polar coordinates (using ellipse 2)

  a = e2.a, b = e2.b;
  var precision = 4 / a;  // results accurate to around 4 metres

  eSq = (a*a - b*b) / (a*a);
  var p = Math.sqrt(x2*x2 + y2*y2);
  var phi = Math.atan2(z2, p*(1-eSq)), phiP = 2*Math.PI;
  while (Math.abs(phi-phiP) > precision) {
    nu = a / Math.sqrt(1 - eSq*Math.sin(phi)*Math.sin(phi));
    phiP = phi;
    phi = Math.atan2(z2 + eSq*nu*Math.sin(phi), p);
  }
  var lambda = Math.atan2(y2, x2);
  H = p/Math.cos(phi) - nu;

  return new LatLon(phi.toDeg(), lambda.toDeg(), H);
}

function getGrid(x1, x2){
        pWGS = new LatLon(x1.parseDeg(),x2.parseDeg());
        pOSGB = convertWGS84toOSGB36(pWGS);
        return LatLongToOSGrid(pOSGB)
}

function OLatLng(lat, lon)
{
   this.lat = lat;
   this.lon = lon;
}

function WGS84LatLng(lat, lon)
{
	this.lat = lat;
	this.lon = lon;
}

function ONorthEast(east, north)
{
	this.north = north;
	this.east = east;
}

function ORect(bottomLeft, topRight)
{
	this.bl = bottomLeft;
	this.tr = topRight;
}

///////////////////////////////////////////////////////////////////////////////
//======================================================================
// based on jscalculators and Bill Chadwick code
//thank you all
// ======================================================================
var deg2rad = Math.PI / 180;
var rad2deg = 180.0 / Math.PI;
var pi = Math.PI;
var locres = 6; //locator characters


//===================================================================
//functions for myform2
//===================================================================
//===================================================================
function transform(lat, lon, a, e, h, a2, e2, xp, yp, zp, xr, yr, zr, s)
{
  // convert to cartesian; lat, lon are radians
  sf = s * 0.000001;
  v = a / (sqrt(1 - (e *(sin(lat) * sin(lat)))));
  x = (v + h) * cos(lat) * cos(lon);
  y = (v + h) * cos(lat) * sin(lon);
  z = ((1 - e) * v + h) * sin(lat);

  xrot = (xr / 3600) * deg2rad;
  yrot = (yr / 3600) * deg2rad;
  zrot = (zr / 3600) * deg2rad;

  hx = x + (x * sf) - (y * zrot) + (z * yrot) + xp;
  hy = (x * zrot) + y + (y * sf) - (z * xrot) + yp;
  hz = (-1 * x * yrot) + (y * xrot) + z + (z * sf) + zp;

  // Convert back to lat, lon
  lon = atan(hy / hx);
  p = sqrt((hx * hx) + (hy * hy));
  lat = atan(hz / (p * (1 - e2)));
  v = a2 / (sqrt(1 - e2 * (sin(lat) * sin(lat))));
  errvalue = 1.0;
  lat0 = 0;
  while (errvalue > 0.001)
  {
    lat0 = atan((hz + e2 * v * sin(lat)) / p);
    errvalue = abs(lat0 - lat);
    lat = lat0;
  }
  h = p / cos(lat) - v;

  return new OLatLng(lat, lon );
}
//===================================================================
function ll2ne(lat, lon, grid)
{
// converts Lat/Lon OSGB to Easting and Northing: inputs to radians.
  var phi = lat * deg2rad; // convert latitude to radians
  var lam = lon * deg2rad; // convert longitude to radians
 // grid = grid;       // British, Irish, Channel

  if (grid == "1")
  {
    a = 6377563.396;       // OSGB semi-major
    b = 6356256.91;        // OSGB semi-minor
    e0 = 400000;           // easting of false origin
    n0 = -100000;          // northing of false origin
    f0 = 0.9996012717;     // OSGB scale factor on central meridian
    e2 = 0.0066705397616;  // OSGB eccentricity squared
    lam0 = -0.034906585039886591;  // OSGB false east
    phi0 = 0.85521133347722145;    // OSGB false north

  }
  if (grid == "2")
  {
    a = 6377340.189;       // OSGBI semi-major
    b = 6356034.447;        // OSGBI semi-minor
    e0 = 200000;           // easting of false origin
    n0 = 250000;          // northing of false origin
    f0 = 1.000035;     // OSGBI scale factor on central meridian
    e2 = 0.00667054015;  // OSGBI eccentricity squared
    lam0 = -0.13962634015954636615389526147909;  // OSGBI false east
    phi0 = 0.93375114981696632365417456114141;    // OSGBI false north

  }
  if (grid == "3")
  {
    a = 6378388.000;       // INT24 ED50 semi-major
    b = 6356911.946;       // INT24 ED50 semi-minor
    e0 = 500000;           // easting of false origin
    n0 = 0;                // northing of false origin
    f0 = 0.9996;           // INT24 ED50 scale factor on central meridian
    e2 = 0.0067226700223333;  // INT24 ED50 eccentricity squared
    lam0 = -0.0523598775598;  // INT24 ED50 false east
    phi0 = 0 * deg2rad;    // INT24 ED50 false north

  }
  try{
  var af0 = a * f0;
  var bf0 = b * f0;

  // easting
  var slat2 = sin(phi) * sin(phi);
  var nu = af0 / (sqrt(1 - (e2 * (slat2))));
  var rho = (nu * (1 - e2)) / (1 - (e2 * slat2));
  var eta2 = (nu / rho) - 1;
  var p = lam - lam0; //LAM-LAM0
  var IV = nu * cos(phi);
  var clat3 = cos(phi) * cos(phi) * cos(phi);
  var tlat2 = tan(phi) * tan(phi);
  var V = (nu / 6) * clat3 * ((nu / rho) - tlat2);
  var clat5 = pow(cos(phi), 5);
  var tlat4 = pow(tan(phi), 4);
  var VI = (nu / 120) * clat5 * ((5 - (18 * tlat2)) + tlat4 + (14 * eta2) - (58 * tlat2 * eta2));
  east = e0 + (p * IV) + (pow(p, 3) * V) + (pow(p, 5) * VI);

  // northing
  var n = (af0 - bf0) / (af0 + bf0);
  var M = Marc(bf0, n, phi0, phi);
  var I = M + (n0);
  var II = (nu / 2) * sin(phi) * cos(phi);
  var III = ((nu / 24) * sin(phi) * clat3) * (5 - tlat2 + (9 * eta2));
  var IIIA = ((nu / 720) * sin(phi) * clat5) * (61 - (58 * tlat2) + tlat4);
  north = I + ((p * p) * II) + (pow(p, 4) * III) + (pow(p, 6) * IIIA);


  // make whole number values
  east = Math.round(east); // round to whole number of meters
  north = Math.round(north);
  }catch (e){
  }
  return new ONorthEast(east, north);
}
//===================================================================
function ne2ll(east, north, grid)  // (east, north, grid)
{
// converts NGR easting and nothing to lat, lon.
// input metres, output radians
  var nX = Number(north);
  var eX = Number(east);

  if (grid == 0)  // no grid specified
    grid = 1;  // default British

  // validate - check all are numbers and something is entered
  if ((String(nX) == "NaN") || (String(eX) == "NaN") || (north.length == 0) || (east.length == 0))
  {
    return;
  }

  if (grid == 2)  // Irish
  {
    a = 6377340.189;       // OSGBI semi-major
    b = 6356034.447;        // OSGBI semi-minor
    e0 = 200000;           // easting of false origin
    n0 = 250000;          // northing of false origin
    f0 = 1.000035;     // OSGBI scale factor on central meridian
    e2 = 0.00667054015;  // OSGBI eccentricity squared
    lam0 = -0.13962634015954636615389526147909;  // OSGBI false east
    phi0 = 0.93375114981696632365417456114141;    // OSGBI false north
  //  en2ngr(round(eX), round(nX), "2");
  }
   if (grid == 1)  // British
  {
    a = 6377563.396;       // OSI semi-major
    b = 6356256.91;        // OSI semi-minor
    e0 = 400000;           // easting of false origin
    n0 = -100000;          // northing of false origin
    f0 = 0.9996012717;     // OSI scale factor on central meridian
    e2 = 0.0066705397616;  // OSI eccentricity squared
    lam0 = -0.034906585039886591;  // OSI false east
    phi0 = 0.85521133347722145;    // OSI false north
    //en2ngr(round(eX), round(nX), "1");
  }
  if (grid == 3)   // Channel Is
  {
    a = 6378388.000;       // INT24 ED50 semi-major
    b = 6356911.946;       // INT24 ED50 semi-minor
    e0 = 500000;           // easting of false origin
    n0 = 0;                // northing of false origin
    f0 = 0.9996;           // INT24 ED50 scale factor on central meridian
    e2 = 0.0067226700223333;  // INT24 ED50 eccentricity squared
    lam0 = -0.0523598775598;  // INT24 ED50 false east
    phi0 = 0 * deg2rad;    // INT24 ED50 false north
    //en2ngr(round(eX), round(nX), "3");
  }
  var af0 = a * f0;
  var bf0 = b * f0;
  var n = (af0 - bf0) / (af0 + bf0);
  var Et = east - e0;
  var phid = InitialLat(north, n0, af0, phi0, n, bf0);
  var nu = af0 / (sqrt(1 - (e2 * (sin(phid) * sin(phid)))));
  var rho = (nu * (1 - e2)) / (1 - (e2 * (sin(phid)) * (sin(phid))));
  var eta2 = (nu / rho) - 1;
  var tlat2 = tan(phid) * tan(phid);
  var tlat4 = pow(tan(phid), 4);
  var tlat6 = pow(tan(phid), 6);
  var clatm1 = pow(cos(phid), -1);
  var VII = tan(phid) / (2 * rho * nu);
  var VIII = (tan(phid) / (24 * rho * (nu * nu * nu))) * (5 + (3 * tlat2) + eta2 - (9 * eta2 * tlat2));
  var IX = ((tan(phid)) / (720 * rho * pow(nu, 5))) * (61 + (90 * tlat2) + (45 * tlat4));
  var phip = (phid - ((Et * Et) * VII) + (pow(Et, 4) * VIII) - (pow(Et, 6) * IX));
  var X = pow(cos(phid), -1) / nu;
  var XI = (clatm1 / (6 * (nu * nu * nu))) * ((nu / rho) + (2 * (tlat2)));
  var XII = (clatm1 / (120 * pow(nu, 5))) * (5 + (28 * tlat2) + (24 * tlat4));
  var XIIA = clatm1 / (5040 * pow(nu, 7)) * (61 + (662 * tlat2) + (1320 * tlat4) + (720 * tlat6));
  var lambdap = (lam0 + (Et * X) - ((Et * Et * Et) * XI) + (pow(Et, 5) * XII) - (pow(Et, 7) * XIIA));

  var geo = convert_to_wgs(grid, phip, lambdap);

  var lat = geo.lat * rad2deg;
  var lon = geo.lon * rad2deg;

  return new OLatLng(lat, lon);
}
//===================================================================
function Marc(bf0, n, phi0, phi)
{
  var Marc = bf0 * (((1 + n + ((5 / 4) * (n * n)) + ((5 / 4) * (n * n * n))) * (phi - phi0))
    - (((3 * n) + (3 * (n * n)) + ((21 / 8) * (n * n * n))) * (sin(phi - phi0)) * (cos(phi + phi0)))
    + ((((15 / 8) * (n * n)) + ((15 / 8) * (n * n * n))) * (sin(2 * (phi - phi0))) * (cos(2 * (phi + phi0))))
    - (((35 / 24) * (n * n * n)) * (sin(3 * (phi - phi0))) * (cos(3 * (phi + phi0)))));
  return(Marc);
}
//===================================================================
function InitialLat(north, n0, af0, phi0, n, bf0)
{
  var phi1 = ((north - n0) / af0) + phi0;
  var M = Marc(bf0, n, phi0, phi1);
  var phi2 = ((north - n0 - M) / af0) + phi1;
  var ind = 0;
  while ((abs(north - n0 - M) > 0.00001) && (ind < 20))  // max 20 iterations in case of error
  {
	ind = ind + 1;
	phi2 = ((north - n0 - M) / af0) + phi1;
    M = Marc(bf0, n, phi0, phi2);
    phi1 = phi2;
  }
  return(phi2);
}
//===================================================================
function choose_where(lat, lon)
{
var where = "";
if ((lat >= 49) && (lat <= 50) && (lon >= -3) && (lon <= -2))
  where = "3"
else if ((lon <= -6) && (lon > -11) && (lat > 51.242) && (lat <= 54))
  where = "2";
else if ((lon < -5.333) && (lon >- 11) && (lat >= 54) && (lat <= 55))
  where = "2";
else if ((lon < -5.9) && (lon >= -9) && (lat >= 55) && (lat < 55.5))
  where = "2";
else if ((lon < 1.8) && (lon > -9) && (lat > 49.8166) && (lat < 61))
  where = "1";
return where;
}




//===================================================================
function convert_to_wgs(grid, phip, lambdap)
{
  var WGS84_AXIS = 6378137;
  var WGS84_ECCENTRIC = 0.00669438037928458;
  var OSGB_AXIS = 6377563.396;
  var OSGB_ECCENTRIC = 0.0066705397616;
  var IRISH_AXIS = 6377340.189;
  var IRISH_ECCENTRIC = 0.00667054015;
  var INT24_AXIS = 6378388.000;
  var INT24_ECCENTRIC = 0.0067226700223333;
  var height = 10;  // dummy height
  var geo = null;
  if (grid == 1)
  {
    geo = transform(phip, lambdap, OSGB_AXIS, OSGB_ECCENTRIC, height, WGS84_AXIS, WGS84_ECCENTRIC, 446.448, -125.157, 542.06, 0.1502, 0.247, 0.8421, 20.4894);
  }
  if (grid == 2)
  {
    geo = transform(phip, lambdap, IRISH_AXIS, IRISH_ECCENTRIC, height, WGS84_AXIS, WGS84_ECCENTRIC, 482.53, -130.596, 564.557, -1.042, -0.214, -0.631, -8.15);
  }
  if (grid == 3)
  {
    geo = transform(phip, lambdap, INT24_AXIS, INT24_ECCENTRIC, height,  WGS84_AXIS, WGS84_ECCENTRIC, -83.901, -98.127, -118.635, 0, 0, 0, 0);
  }
  return(geo);
}
function convert_to_ogb(grid, phip, lambdap)
{
  var WGS84_AXIS = 6378137;
  var WGS84_ECCENTRIC = 0.00669438037928458;
  var OSGB_AXIS = 6377563.396;
  var OSGB_ECCENTRIC = 0.0066705397616;
  var IRISH_AXIS = 6377340.189;
  var IRISH_ECCENTRIC = 0.00667054015;
  var INT24_AXIS = 6378388.000;
  var INT24_ECCENTRIC = 0.0067226700223333;
  var height = 0;  // dummy height
  var geo = null;
  if (grid == 1)
  {
    geo = transform(phip, lambdap, WGS84_AXIS, WGS84_ECCENTRIC , height, OSGB_AXIS, OSGB_ECCENTRIC, -446.448, 125.157, -542.06, -0.1502, -0.247, -0.8421, 20.4894);
  }
  if (grid == 2)
  {
    geo = transform(phip, lambdap, WGS84_AXIS, WGS84_ECCENTRIC , height, IRISH_AXIS, IRISH_ECCENTRIC, -482.53, 130.596, -564.557, 1.042, 0.214, 0.631, -8.15);
  }
  if (grid == 3)
  {
    geo = transform(phip, lambdap, WGS84_AXIS, WGS84_ECCENTRIC , height,  INT24_AXIS, INT24_ECCENTRIC, 83.901, 98.127, 118.635, 0, 0, 0, 0);
  }
  return(geo);
}

//===================================================================
function en2ngr(east, north, grid)
{
  var nstr = String(north);
  var estr = String(east);
  while (nstr.length < 6) {nstr = "0" + nstr};
  while (estr.length < 6) {estr = "0" + estr};

  if (grid == "3")
  {

    if (north > 5500000)
    {
      var ngr = "WA " + estr.substring(1, 7) + " " + nstr.substring(2, 7);

    }
    else if (north < 5500000)
    {
      var ngr = "WV " + estr.substring(1, 7) + " " + nstr.substring(2, 7);

    }
  }

  if ((grid == "1") || (grid == "2"))
  {
    var eX = east / 500000;
    var nX = north / 500000;
    var tmp = floor(eX)-5.0 * floor(nX)+17.0;
    nX = 5 * (nX - floor(nX));
    eX = 20 - 5.0 * floor(nX) + floor(5.0 * (eX - floor(eX)));
    if (eX > 7.5)
      eX = eX + 1;
    if (tmp > 7.5)
      tmp = tmp + 1;
    var eing = estr;
    var ning = nstr;

    var l = eing.length;
    eing = eing.substring(l - 5, l);
    l = ning.length;
    ning = ning.substring(l - 5, l);
    ngr = chr(tmp + 65) + chr(eX + 65) + " " + eing + " " + ning;
  }


  if (grid == "2")
  {
    l = ngr.length;
    ngr = "I"+ngr.substring(1, l);
  }
  return ngr;
}

//===================================================================
function conv_ngr_to_ings(ngr, grid)
{
  var ngr = ngr.toUpperCase();
  var lenngr = ngr.length;

  var i = 0;
  var t = "";

  // replace spaces
  while(i < lenngr)
  {
    if (ngr.charAt(i) != " ")
      t = t + ngr.charAt(i);
  i++;
  }

  // check length, then add leading space if odd - Irish?
  lenngr = t.length;
  if (lenngr >12)
  {
    return;
  }

  var temp = round(lenngr/2)*2;
  if (lenngr != temp)
    ngr = " " + t;  // possible Irish NGR
  else
    ngr = t;

  // now format to 12 digits
  if (ngr.length == 4)
    t = ngr.substring(0,3) + "5000" + ngr.substring(3) + "5000";
  if (ngr.length == 6)
    t = ngr.substring(0,4) + "000" + ngr.substring(4) + "000";
  if (ngr.length == 8)
    t = ngr.substring(0,5) + "00" + ngr.substring(5) + "00";
  if (ngr.length == 10)
    t = ngr.substring(0,6) + "0" + ngr.substring(6) + "0";
  if (ngr.length == 12)
    t = ngr;
  ngr = t;

  // check ngr has two letters and 10 numbers
  var myRegExp = /[ A-Z]{2}[0-9]{10}/;
  if (! myRegExp.test(ngr))
  {
    return;
  }

  if (grid == "3")
  {
    var eci = "5" + ngr.substring(2,7);
    if (ngr.charAt(1) == "V")
      var nci = "54" + ngr.substring(7,12);
    else
      var nci = "55" + ngr.substring(7,12);
    var east = Number(eci);
    var north = Number(nci);
  }
  if ((grid == "1") || (grid == "2"))
  {
    var irish = 0;  //not Irish
    var north = Number(ngr.substring(7));
    var east = Number(ngr.substring(2,7));

    var t1 = ngr.charCodeAt(0) - 65;
    if (t1 < 0)
    {
      t1 = 18;    // S for Irish 83-65
      irish = 1;
    }

    if (t1 > 8)
      t1 = t1 -1;
    var t2 = floor(t1 / 5);
    north = north + 500000 * (3 - t2);
    east = east + 500000 * (t1 - 5 * t2 - 2);

    t1 = ngr.charCodeAt(1) - 65;
    if (t1 > 8)
      t1 = t1 - 1;
    t2 = floor(t1 / 5);
    north = north + 100000 * ( 4 - t2);
    east = east + 100000 * ( t1 - 5 * t2);
  }
  return new ONorthEast(east, north);

}

function ngr2wgs(ngr, grid){
  var ings = conv_ngr_to_ings(ngr, grid);
  var ll = ne2ll(ings.east, ings.north, grid);

  return ll;

}

//
//convert WGS84 Latitude and Longitude to Ordnance Survey 1936 Latitude Longitude

function wgs2ogb(WGlat, WGlon, grid)
{
  //first off convert to radians
  var radWGlat = WGlat * deg2rad;
  var radWGlon = WGlon * deg2rad;

  var geo = convert_to_ogb(grid, radWGlat, radWGlon);
  //convert back to degrees
  var newLat = geo.lat * rad2deg;
  var newLon = geo.lon * rad2deg;

  return new OLatLng(newLat, newLon);

}

//return rect in OSGB N/E coords that encloses the given wgs84 rect
function enclosingOsgbRect(WGleft,WGbottom,WGtop,WGright, grid){

	var blOGB = wgs2ogb(WGbottom,WGleft,grid);
	var trOGB = wgs2ogb(WGtop,WGright,grid);
	var brOGB = wgs2ogb(WGbottom,WGright,grid);
	var tlOGB = wgs2ogb(WGtop,WGleft,grid);

	var blEN = ll2ne(blOGB.lat,blOGB.lon,grid);
	var trEN = ll2ne(trOGB.lat,trOGB.lon,grid);
	var brEN = ll2ne(brOGB.lat,brOGB.lon,grid);
	var tlEN = ll2ne(tlOGB.lat,tlOGB.lon,grid);

	var e = Math.min(blEN.east,tlEN.east);
	var w = Math.max(brEN.east,trEN.east);
	var s = Math.min(blEN.north,brEN.north);
	var n = Math.max(trEN.north,tlEN.north);

	return new ORect(new ONorthEast(e,s),new ONorthEast(w,n));
}


//===================================================================


// Maths functions

function mod(y, x)
{
if (y >= 0)
  return y - x * floor(y / x);
else
  return y + x * (floor(-y / x) + 1.0);
}

function atan2(y, x)
{
	return Math.atan2(y, x);
}

function sqrt(x)
{
	return Math.sqrt(x);
}

function tan(x)
{
	return Math.tan(x);
}

function sin(x)
{
	return Math.sin(x);
}

function cos(x)
{
	return Math.cos(x);
}

function acos(x)
{
	return Math.acos(x);
}

function floor(x)
{
	return Math.floor(x);
}

function round(x)
{
	return Math.round(x);
}

function ceil(x)
{
	return Math.ceil(x)
}

function ln(x)
{
	return Math.log(x);
}

function abs(x)
{
	return Math.abs(x);
}

function pow(x, y)
{
	return Math.pow(x, y);
}

function atan(x)
{
	return Math.atan(x);
}

function chr(x)
{
return String.fromCharCode(x);
}

function round(x)
{
	return Math.round(x);
}

