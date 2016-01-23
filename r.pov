#default{ finish{ ambient 0.1 diffuse 0.9 conserve_energy}}
global_settings { max_trace_level 10 }
background { color rgb <1, 1, 1> }

//#declare cam_pos = <20, -10, -8>;
//camera { location cam_pos look_at  <0,-6.2,0> angle 5}
#declare cam_pos = <55, 10, -35>;
camera { location cam_pos look_at  <0,0,0> angle 22}

light_source{<500,2500,0> color 0.9}   // sun light
light_source{ cam_pos  color rgb<0.9,0.9,1>*0.8}  // flash light

/******************************************************/
#macro copper()
    material{
      texture {
        pigment {rgb<1,0.5,0.3>}
        normal { bumps 0.1 scale 0.4}
        finish { diffuse 0.9 reflection 0.1
                 specular 0.4 roughness 0.3 phong 1 phong_size 40}
      }
    }
#end
#macro copper1()
    material{
      texture {
        pigment {rgb<0.7,0.5,0.3>}
        normal { bumps 0.2 scale 0.2}
        finish { diffuse 0.9 reflection 0.1
                 specular 0.4 roughness 0.3 phong 1 phong_size 40}
      }
    }
#end

#macro copper2()
    material{
      texture {
        pigment {rgb<0.5,0.5,0.5>}
        normal { bumps 0.05 scale 0.2}
        finish { diffuse 0.9 reflection 0.1
                 specular 0.4 roughness 0.3 phong 1 phong_size 40}
      }
    }
#end

#macro glass()
    material{
      texture {
        pigment{ rgbf <0.98, 0.98, 0.98, 0.9> }
        normal { bumps 0.05 scale 0.09}
        finish { diffuse 0.1 reflection 0.2
                 specular 0.8 roughness 0.0003 phong 1 phong_size 400}
      }
      interior{ ior 1.5 caustics 0.5}
    }
#end

/*************************************************/

/* The main ring*/
#local RING_RAD = 5;
#local RING_TH = 0.08;
#local RING_GAP = 1.2;
#local RING_GAP1 = 2.0;
#local RING_ANG = (180-104)/2;

/* The central ring*/
#local CENT_RAD = 1;

/* the arrow */
#local ARROW_LEN = 7;    // short arrow length
#local ARROW_RAD = 0.25;  // arrow radius
#local ARROW_HRD = 0.5;  // arrow head radius
#local ARROW_HLN = 1.0;  // arrow head length

#local HOLE_RAD =0.12;  // hole radius

#local ARROW_GAP = 0.9; // gap between arrows
#local ARROW_ANG = 40;  // angle between arrows
#local SMRING = 2.5*RING_TH; // small ring rad

#local ARM_PT  = RING_RAD/3;  // point where arms are connected
#local ARM_LEN = RING_RAD/2; // length of arms


/*************************************************/

#macro main_ring(RAD)
merge{
  #local BBOX = RAD+1;
  difference{
    torus{RAD,RING_TH}
    box{<BBOX,-1,0>,<-BBOX,1,BBOX>}
    rotate (180-104)/2*x
    translate -RING_GAP/2*z
  }
  difference{
    torus{RAD,RING_TH}
    box{<BBOX,-1,0>,<-BBOX,1,-BBOX>}
    rotate -RING_ANG*x
    translate +RING_GAP/2*z
  }
  merge{
    cylinder{<0,0,-1>, <0,0,1>, RING_TH translate RAD*x}
    cylinder{<0,0,-1>, <0,0,1>, RING_TH translate -RAD*x}
    scale RING_GAP/2*z
  }
/*  merge{
    #local A = pi/180*RING_ANG;
    #local H = (RING_GAP1-RING_GAP)/2*tan(A);
    #local L = (RING_GAP1-RING_GAP)/2/sin(A);
    #local D = sqrt(RAD*RAD-L*L);
    cylinder{<0,0,-1>, <0,0,1>, RING_TH translate D*x}
    cylinder{<0,0,-1>, <0,0,1>, RING_TH translate -D*x}
    scale RING_GAP1/2*z
    translate H*y
  }
*/
 torus{ARROW_RAD+RING_TH,RING_TH
       rotate 90*z
       translate (RAD-RING_TH)*x}
 torus{ARROW_RAD+RING_TH,RING_TH
       rotate 90*z
       translate -(RAD-RING_TH)*x}
}
#end

/*************************************************/
#macro arm(S)
  merge{
    #local LL1 = ARROW_RAD+1.5*RING_TH;
    #local LL2 = 1.5*RING_TH;
    cylinder{<SL1,SL2,S*LL1>,<SLPOS,0,S*LL2>,RING_TH} // arm itself
    cylinder{<SL1,SL2,S*LL1>,<SL1,SL2,-S*LL1>,RING_TH} // bending in the arrow
    cylinder{<SLPOS,0,S*LL2>,<SLPOS,0,-S*LL2>,RING_TH}
    sphere{<SL1,SL2,-S*LL1>,2*RING_TH}
    sphere{<SLPOS,0,-S*LL2>,2*RING_TH}
    sphere{<SL1,SL2,S*LL1>,RING_TH}
    sphere{<SLPOS,0,S*LL2>,RING_TH}
    rotate S*(90-RING_ANG)*x
    translate S*RING_GAP/2*z
  }
#end

#macro slider()
  #local SL1 = ARM_PT*cos(ARROW_ANG/180*pi);
  #local SL2 = ARM_PT*sin(ARROW_ANG/180*pi);
  #local SLPOS = SL1 + sqrt(ARM_LEN*ARM_LEN - SL2*SL2); // slider position
  merge{
    torus{ARROW_RAD+RING_TH,RING_TH
          rotate 90*z translate SLPOS*x}
    torus{ARROW_RAD+RING_TH,RING_TH
          rotate 90*z translate (SLPOS-2*RING_TH)*x}
    torus{ARROW_RAD+RING_TH,RING_TH
          rotate 90*z translate (SLPOS+2*RING_TH)*x}
    torus{SMRING,RING_TH
          rotate -RING_ANG*x
          translate SLPOS*x+RING_GAP/2*z}
    torus{SMRING,RING_TH
          rotate RING_ANG*x
          translate SLPOS*x-RING_GAP/2*z }
    // arms
    arm(+1)
    arm(-1)
  }
#end

#macro small_ring()
  torus{  ARROW_HLN,RING_TH
          translate z*(ARROW_LEN+ARROW_HLN/2)
          rotate (90-ARROW_ANG)*y
          rotate -RING_ANG*x
          translate z*RING_GAP/2
  }
#end

/*************************************************/

#macro arrow(L)
  #local RR = 0.1;  // rounding radius, mm
  merge{
    cylinder{<0,0,0>,<0,L-ARROW_HLN+RR,0>,ARROW_RAD}
    sphere{<0,0,0>, ARROW_RAD}
    sphere{<0,L,0>,RR}
    cone{<0,L-ARROW_HLN,0>, ARROW_HRD-RR <0,L,0>,0}
    cone{<0,L-ARROW_HLN+RR,0>,ARROW_HRD,<0,L,0>,RR}
    torus{ARROW_HRD-RR,RR translate y*(L-ARROW_HLN+RR)}
    difference{
      cylinder{<0,L-ARROW_HLN-RR,0>,<0,L-ARROW_HLN+RR,0>, ARROW_RAD+RR}
      torus{ARROW_RAD+RR,RR translate y*(L-ARROW_HLN-RR)}
    }
  }
#end
#macro vectorN()
  difference{
    object{arrow(2*ARROW_LEN) translate -ARROW_LEN*y rotate -90*z}
    cylinder{<0,0,-ARROW_HRD>, <0,0,ARROW_HRD>, HOLE_RAD translate -RING_RAD*x}
    cylinder{<0,0,-ARROW_HRD>, <0,0,ARROW_HRD>, HOLE_RAD translate +RING_RAD*x}
  }
#end
#macro vectorLS()
  difference{
    object{arrow(ARROW_LEN-ARROW_GAP) translate ARROW_GAP*y}
    cylinder{<-ARROW_HRD,0,0>, <ARROW_HRD,0,0>, HOLE_RAD translate RING_RAD*y}
    cylinder{<-ARROW_HRD,0,0>, <ARROW_HRD,0,0>, HOLE_RAD translate CENT_RAD*y}
    cylinder{<0,0,-ARROW_HRD>, <0,0,ARROW_HRD>, HOLE_RAD translate ARM_PT*y}
    cylinder{<-ARROW_HRD,0,0>, <ARROW_HRD,0,0>, HOLE_RAD translate (ARROW_LEN+ARROW_HLN/2)*y}
    rotate -(90-ARROW_ANG)*z
  }
#end
#macro vectorL()
  object{
    vectorLS()
    rotate (90-RING_ANG)*x 
    translate +RING_GAP/2*z
  }
#end
#macro vectorS()
  object{
    vectorLS()
    rotate -(90-RING_ANG)*x
    translate -RING_GAP/2*z
  }
#end


//arrow()
//torus{0.8,0.05 rotate 90*x translate -5.3*y copper1()}

merge{
merge{
  main_ring(RING_RAD)
  main_ring(CENT_RAD)
  slider()
  small_ring()
  copper1()
}
merge{
  vectorN()
  vectorL()
  vectorS()
//  pigment{rgb<1,0,0>}
  glass()
//  copper2()
}
  rotate -(90-RING_ANG)*x
  rotate (90-ARROW_ANG)*z
  rotate -55*y
}