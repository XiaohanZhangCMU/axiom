#include "linalg3.h"


Vector3 operator + (const Vector3 &a, const Vector3 &b){ 
 return Vector3(a.x+b.x, a.y+b.y, a.z+b.z);
}

Vector3 operator - (const Vector3 &a, const Vector3 &b){ 
 return Vector3(a.x-b.x, a.y-b.y, a.z-b.z);
}

Vector3 operator * (const Vector3 &dv, const double d) {
 return Vector3(dv.x*d, dv.y*d, dv.z*d);
}

Vector3 operator * (const double d, const Vector3 &dv) {
 return Vector3(dv.x*d, dv.y*d, dv.z*d);
}

Vector3 operator / (const Vector3 &dv, const double d) {
        return Vector3(dv.x/d, dv.y/d, dv.z/d);
    }


