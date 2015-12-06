
var minX = -10,
    minY = -5,
    maxX = -1e-2,
    maxY = 5,
    dx   = 0.13,
    dt   = 0.05,
    h  = 0.05,

    l = 2,
    sigma = 4,
    f = 0.4,
    a = 7,
    k = 1,

    exp = math.exp,
    sin = math.sin,
    cos = math.cos,
    add = math.add,
    mult = math.multiply,
    cos = math.cos,
    i = math.complex(0,1),
    re = math.re
    ;

var ss = new SimShim(
  document.getElementById("plot"),
  {
    lightIntensity: 1.0
  }
);

/////////////////////////////////////////////////////////////////////// ETA

function eta (x,y,t) {
  // return exp( add( mult(i,l*y - sigma*t), f*l*x/sigma ) ).re;
  return cos( l*y - sigma*t )*exp( f*l*x/sigma );
}

function createMesh(minx,miny,maxx,maxy,step,fn) {
  // this creates a 2D array of z values
  var sb = [];
  for (var i = minx; i < maxx; i+=step) {
    var row = [];
    for (var j = miny; j < maxy; j+=step) {
      row.push( fn(i,j) );
    }
    sb.push( row );
  }
  return sb;
}

function iterMesh(mesh, t) {
  // return relaxation iteration and max value of V'ij - Vij
  var maxVal = 0,
      hgt    = mesh.length,
      wdt    = mesh[0].length,
      m      = [];

  for (var j = 0; j < hgt; j++) {
    var row = [];
    for (var i = 0; i < wdt; i++) {
      var x = minX + (i/wdt)*(maxX - minX),
          y = minY + (j/hgt)*(maxY - minY),
          m_ji = eta(x,y,t);
      row.push(m_ji);
    }
    m.push(row);
  }
  return m;
}

var eta_plt = {
  "type"     : "surfaceplot",
  "animated" : true,
  "minX"     : minX,
  "minY"     : minY,
  "maxX"     : maxX,
  "maxY"     : maxY,
  "start"    : 0,
  "mesh"     : createMesh(minX, minY, maxX, maxY, dx, eta),
  "t"        : 0,
  "next"     : function () {
    this.t += dt;
    return iterMesh(this.mesh, this.t);
  }
};

ss.addPlot( eta_plt, {"color": "#00c0ff"} );

//////////////////////////////////////////////////////////////////// PATHLINES

/*\
|*|  Runge Kutta
\*/

function rk4 (xyz, f, t) {
  var k1,k2,k3,k4,xyz1,xyz2,xyz3,xyz4;
  // RK method
  k1   = sm(   h,   f( xyz, t ) );
  xyz1 = add3( xyz, sm( 0.5, k1 ));
  k2   = sm(   h,   f( xyz1, t + h/2));
  xyz2 = add3( xyz, sm( 0.5, k2 ));
  k3   = sm(   h,   f( xyz2, t + h/2 ));
  xyz3 = add3( xyz, k3 );
  k4   = sm(   h,   f( xyz3, t+h ));
  // Compute estimated component
  xyz4 = add3(
      xyz,
      sm( 1/6, k1 ),
      sm( 1/3, k2 ),
      sm( 1/3, k3 ),
      sm( 1/6, k4 )
  );
  return xyz4;
}

function sm (c, arr) {
    // scalar multiplication
    return arr.map(function (x) { return c*x });
}

function add3 (/* variable number of arrays */) {
    // add 3 dimensional arrays
    var arrays = arguments;
    var arr = [];
    for (var i = 0; i < 3; i++) {
        arr_i = 0;
        for (var j = 0; j < arrays.length; j++) {
            arr_i += arrays[j][i];
        }
        arr.push(arr_i);
    };
    return arr;
}

function kelvinDeep(xyz,t) {
  var x = xyz[0],
      y = xyz[1],
      z = xyz[2],
      foo = a*sigma*exp(k*z);
  return [
    0,
    foo*cos(k*x-sigma*t),
    foo*sin(k*x-sigma*t)
  ];
}

/*\
|*|  ADD PARTICLES
\*/

function addParticle() {
  var x = -1 - 8*Math.random(),
      y = -3 + 5*Math.random(),
      z = -2.8 - Math.random();
      plotID = ss.addPlot({
    // required properties (used by SimShim)
    "type": "lineplot",
    "animated": true,
    "lineLength": 100,
    "xyz": [x,y,z], // initial condition
    // The "next" function must return the next point each frame.
    "next": function () {
      this.posn = rk4(this.posn, kelvinDeep, this.t);
      this.t += dt;
      return this.posn;
    },
    // custom properties
    "posn": [x,y,z],
    "t": 0
  });
  return plotID;
}

plotIDs = [];

for (var i = 0; i < 6; i++) {
  plotIDs.push( addParticle() );
}

// remove out-of-bounds plots
function cull() {
  for (var i = 0; i < plotIDs.length; i++) {
    var pid = plotIDs[i],
        plt = ss.getPlot(pid);
    if (plt.obj.xyz.y > 5 || plt.obj.t > 7) {
      ss.removeById( pid );
      plotIDs.splice(i,1);
      plotIDs.push( addParticle() );
    }
  }
}

setInterval(cull, 1000);

/////////////////////////////////////////////////////////////////////// Custom

var dl = new THREE.DirectionalLight(0xffffff, 0.7);
dl.position.set(-5, 0, 10);
dl.lookAt(new THREE.Vector3(-5,0,0));
ss.plotCtx.scene.add( dl );

//////////////////////////////////////////////////////////////////////// Start

ss.start();
