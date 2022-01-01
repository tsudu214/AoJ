#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <climits> 
#include <vector>
#include <deque>
#include <list>
#include <stack>
#include <queue>
#include <set>
#include <map>
#include <bitset>
#include <string>
#include <cstring>
#include <cmath>
#include <array>
#include <numeric>
#include <unordered_map>
#include <cassert>
#include <random>
#include <chrono>

using namespace std;
using namespace std::chrono;

using ll = long long;
const ll INF = LLONG_MAX;

//------------------------------------------------------------------
// [usage example]
// Clock clock
// ...
// cerr << "proc A [" << clock.Duration<microseconds>().count()/1000. << " msec]" << "\n"; 
class Clock
{
public:
    Clock() {
        Reset();
    }

    void Reset() {
        begin = high_resolution_clock::now();
    }

    template <typename T>
    T Duration() {
        T d = duration_cast<T>(high_resolution_clock::now() - begin);
        begin = high_resolution_clock::now();
        return d;
    }

private:
    high_resolution_clock::time_point begin;
};

//------------------------------------------------------------------
// 2D geometry class
struct Point 
{
    Point() = default;
    Point(double x0, double y0) : x(x0), y(y0) {}

    Point operator+(const Point& q) const { return {x + q.x, y + q.y}; }
    Point operator-(const Point& q) const { return {x - q.x, y - q.y}; }
    Point operator*(double a) const { return {a * x, a * y}; }

    Point& operator+=(const Point& q) { x += q.x; y += q.y; return *this; }
    Point& operator-=(const Point& q) { x -= q.x; y -= q.y; return *this; }
    Point& operator*=(double a) { x *= a; y *= a; return *this; }

    double sqlen() const {
        return x * x + y * y;
    }

    double len() const {
        return sqrt(sqlen());
    }

    bool equal(const Point &q, double eps = 1e-10) const {
        return fabs(x - q.x) < eps && fabs(y - q.y) < eps;
    }

    Point norm() const {
        double a = 1./len();
        return {a*x, a*y};
    }

    double x;
    double y;
};

bool compare_x(const Point &a, const Point &b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

bool compare_y(const Point &a, const Point &b) {
    return a.y < b.y || (a.y == b.y && a.x < b.x);
}

typedef Point Vector;

using Polygon = vector<Point>;

double clamp(double x, double min, double max) 
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

double dot(const Vector& p, const Vector& q) 
{
    return p.x * q.x + p.y * q.y;
}

// outer-product of 2D vectors is scalar
double cross(const Vector& p, const Vector& q) {
    return p.x * q.y - p.y * q.x;
} 

struct Line 
{
    Line(const Point& p0, const Vector& v0) : p(p0), v(v0) {}
    Point p;
    Vector v;  // should be normalized
};

struct Segment 
{
    Segment(const Point& a0, const Point& b0): a(a0), b(b0) {}
    Point a;
    Point b;
    Line line(double* len = nullptr) const {
        if(len) *len = (b - a).len();
        return Line(a, (b - a).norm());
    }
};

struct AABB {
    AABB() : min({HUGE_VAL, HUGE_VAL}), max({-HUGE_VAL, -HUGE_VAL}) {}
    AABB(const Point& m, const Point& M) : min(m), max(M) {}
    AABB(const Polygon& poly) : min({HUGE_VAL, HUGE_VAL}), max({-HUGE_VAL, -HUGE_VAL})
    {
        for (const auto& v : poly) {
            min.x = std::min(min.x, v.x);
            min.y = std::min(min.y, v.y);
            max.x = std::max(max.x, v.x);
            max.y = std::max(max.y, v.y);
        }
    }
    Point min;
    Point max;
};

double param(const Line& l, const Point& p)
{
    return dot(p - l.p, l.v);
}

Point getpoint(const Line& l, double t) 
{
    return l.p + l.v * t;
}

Point project(const Line& l, const Point& p, double* t = nullptr, double* dist = nullptr) 
{
    double th = param(l, p);
    Point h = getpoint(l, th);
    if (dist) *dist = (p - h).len();
    if (t) *t = th;
    return h;
}

Point symmetry(const Line& l, const Point& p)
{
    return p + (project(l, p) - p) * 2;
}

double signedArea(const Point& p0, const Point& p1, const Point& p2)
{
    return cross((p1 - p0), (p2 - p0));
}

// 1:ccw, -1:cw, 0:on-line
int ccw(const Point& p0, const Point& p1, const Point& p2, double eps, double* area = nullptr)
{
    double a = signedArea(p0, p1, p2);
    if (area) *area = a;
    if (a > eps) 
        return 1;
    else if (a < -eps)
        return -1;
    else
        return 0;
}

bool intersect(const Segment& s1, const Segment& s2, double eps, Point* x = nullptr)
{
    const Point& p1 = s1.a;
    const Point& p2 = s1.b;
    const Point& p3 = s2.a;
    const Point& p4 = s2.b;

    double a1, a2, a3, a4;
    int ccw1 = ccw(p1, p2, p3, eps, &a1);
    int ccw2 = ccw(p1, p2, p4, eps, &a2);
    int ccw3 = ccw(p3, p4, p1, eps, &a3);
    int ccw4 = ccw(p3, p4, p2, eps, &a4);

    if (ccw1*ccw2 <= 0 && ccw3*ccw4 <= 0) {
        if (ccw1 == 0 && ccw2 == 0 && ccw3 == 0 && ccw4 == 0) {
            double s = 0;
            double t = (p2 - p1).len();
            Line l(p1, (p2 - p1).norm());
            double m3, m4;
            Point h3 = project(l, p3, &m3);
            Point h4 = project(l, p4, &m4);
            if (m3 > m4) swap(m3, m4);
            if (s < m4 + eps && m3 < t + eps) {
                if (x) {
                    double m = 0.5*(max(s, m3) + min(t, m4));
                    *x = getpoint(l, m);
                }
                return true;
            }
            return false;
        }
        else {
            if (x) {
                if (ccw1 != 0 || ccw2 != 0) {
                    double m = a2 / (a2 - a1);
                    *x = p3 * m + p4 * (1-m);
                }
                else { // ccw3 != 0 || ccw4 != 0
                    double m = a4 / (a4 - a3);
                    *x = p1 * m + p2 * (1-m);
                }
            }
            return true;
        }
    }
    else {
        return false;
    }
}

Point project(const Segment& s, const Point& p, double* t = nullptr, double* dist = nullptr) 
{
    Vector ap = p - s.a;
    Vector ab = s.b - s.a;
    double e = dot(ap, ab);
    double f = dot(ab, ab);

    if (e <= 0) {
        if (t) *t = 0;
        if (dist) *dist = ap.len();
        return s.a;
    }
    if (e >= f) {
        if (t) *t = ab.len();
        if (dist) *dist = (p - s.b).len();
        return s.b;
    }
    return project(s.line(), p, t, dist);
}

bool intersect(const Line& l1, const Line& l2, double eps, double* sx = nullptr, double* tx = nullptr, Point* x = nullptr)
{
    double crs_v = cross(l1.v, l2.v);
    if (fabs(crs_v) < eps) { // parallel
        double dist;
        double t;
        Point h = project(l2, l1.p, &t, &dist);
        if (dist < eps) {
            if (sx) *sx = 0.;
            if (tx) *tx = t;
            if (x) *x = l1.p;
            return true;
        }
        return false;
    }
    double crs_v2 = cross(l2.p - l1.p, l2.v);
    double s = crs_v2 / crs_v;
    double crs_v1 = cross(l2.p - l1.p, l1.v);
    double t = crs_v1 / crs_v;

    Point x1 = getpoint(l1, s);
    Point x2 = getpoint(l2, t);
    assert(x1.equal(x2, eps));

    if (sx) *sx = s;
    if (tx) *tx = t;
    if (x) *x = getpoint(l1, s);

    return true;
}

void closest(const Segment& seg1, const Segment& seg2, double eps, 
    double* dist, 
    pair<double, Point>* sx = nullptr, pair<double, Point>* tx = nullptr)
{
    double smax, tmax;
    Line l1 = seg1.line(&smax);
    Line l2 = seg2.line(&tmax);

    double s, t;
    Point x;
    bool is_intersect = intersect(l1, l2, eps, &s, &t, &x);
    if (is_intersect) {
        double s1 = clamp(s, 0, smax);
        double t1 = clamp(t, 0, tmax);
        if (s1 == s && t1 == t) {
            *dist = 0;
            if (sx) *sx = make_pair(s, x);
            if (tx) *tx = make_pair(t, x);
        }
        else if (s1 == s) {
            Point p2 = getpoint(l2, t1);
            double s2;
            Point p1 = project(seg1, p2, &s2, dist);
            if (sx) *sx = make_pair(s2, p1);
            if (tx) *tx = make_pair(t1, p2);
        }
        else if (t1 == t) {
            Point p1 = getpoint(l1, s1);
            double t2;
            Point p2 = project(seg2, p1, &t2, dist);
            if (sx) *sx = make_pair(s1, p1);
            if (tx) *tx = make_pair(t2, p2);
        } else {
            Point p1 = getpoint(l1, s1);
            double s2, t2;
            Point p2 = project(seg2, p1, &t2, dist);
            p1 = project(seg1, p2, &s2, dist);
            if (sx) *sx = make_pair(s2, p1);
            if (tx) *tx = make_pair(t2, p2);
        }
    }
    else {
        double s1, s2, d1, d2;
        Point h1 = project(seg1, seg2.a, &s1, &d1);
        Point h2 = project(seg1, seg2.b, &s2, &d2);
        if (s1 > s2) swap(s1, s2);
        if (s1 < smax + eps && s2 > -eps) { // overlapping
            double m = 0.5*(max(0., s1) + min(smax, s2));
            Point p1 = getpoint(l1, m);
            double t2;
            Point p2 = project(seg2, p1, &t2, dist);
            if (sx) *sx = make_pair(m, p1);
            if (tx) *tx = make_pair(t2, p2);
        }
        else {
            if (d1 < d2) {
               s = param(seg1.line(), h1);
               *dist = (seg2.a - h1).len();
               if (sx) *sx = make_pair(s, h1);
               if (tx) *tx = make_pair(0.0, seg2.a);
            }
            else {
               s = param(seg1.line(), h2);
               *dist = (seg2.b - h2).len();
               if (sx) *sx = make_pair(s, h2);
               if (tx) *tx = make_pair(tmax, seg2.b);
            }
        }
    }
}

bool is_convex(const Polygon& poly, double eps)
{
    int n = (int)poly.size();
    for (int i = 0; i < n; i++) {
        int prev = (i + n -1) % n;
        int next = (i + 1) % n;
        if (ccw(poly[prev], poly[i], poly[next], eps) < 0) {
            return false;
        }
    }
    return true;
}

static double get_widest_angle_diff( const vector<double>& v)
{
    double ans = 0;
    double w_max = 0;
    int n = (int)v.size();
    for (int i = 0; i < n-1; i++) {
        double w = v[i+1] - v[i];
        if (w > w_max) {
            ans = 0.5 * (v[i+1] + v[i]);
            w_max = w;
        }
    }
    return ans;
}

// 2:inside,  1:on-edge, 0:outside
int contain(const Polygon& poly, const Point& p, double eps)
{
    AABB mmb(poly);
    double diag = (mmb.max - mmb.min).len();

    int n = (int)poly.size();
    for (int i = 0; i < n; i++) {
        int next = (i + 1) % n;
        double t, dist;
        Point h = project(Segment(poly[i], poly[next]), p, &t, &dist);
        if (dist < eps) {
            return 1;
        }
    }

    vector<double> arg; // angle
    for (const auto& q : poly) {
        Vector v = q - p;
        double a = atan2(v.y, v.x);
        arg.push_back(a);
    }
    sort(arg.begin(), arg.end());

    double La = get_widest_angle_diff(arg);
    Vector vv(cos(La), sin(La));
    Point q = p + vv * diag;
    Segment ray(p, q);

    int count = 0;
    for (int i = 0; i < n; i++) {
        int next = (i + 1) % n;
        Segment e(poly[i], poly[next]);

        if (intersect(ray, e, eps)) {
            count++;
        }
    }

    if (count % 2 == 0) {
        return 0;
    } else {
        return 2;
    }
}

double area(const Polygon& poly)
{
    int n = (int)poly.size();
    double s = 0;
    for (int i = 1; i < n-1; i++) {
        s += 0.5 * signedArea(poly[0], poly[i], poly[(i+1)%n]);
    }
    return s;
}

Polygon convex_hull(const Polygon& poly, double eps)
{
    vector<Point> s = poly;
    int n = (int)s.size();
    if (n < 3) return s;

    sort(s.begin(), s.end(), compare_y);

    vector<Point> up, down;
    up.reserve(n);
    down.reserve(n);
    up.push_back(s[0]);
    up.push_back(s[1]);
    for (int i = 2; i < n; i++) {
        for (int j = (int)up.size()-1; j >= 1; j--) {
            if (ccw(up[j-1], up[j], s[i], eps) == 1) {
                up.pop_back();
            } else {
                break;
            }
        }
        up.push_back(s[i]);
    }

    down.push_back(s[n-1]);
    down.push_back(s[n-2]);
    for (int i = n-3; i >= 0; i--) {
        for (int j = (int)down.size()-1; j >= 1; j--) {
            if (ccw(down[j-1], down[j], s[i], eps) == 1) {
                down.pop_back();
            } else {
                break;
            }
        }
        down.push_back(s[i]);
    }

    reverse(down.begin(), down.end());
    for (int i = (int)up.size()-2; i >= 1; i--) {
        down.push_back(up[i]);
    }

    return down;
}

double convex_diameter(const Polygon& poly)
{
    int imin = 0, imax = 0;
    int n = (int)poly.size();
    double ymin = HUGE_VAL, ymax = -HUGE_VAL;
    for (int i = 0; i < n; i++) {
        if (poly[i].y < ymin) { 
            ymin = poly[i].y;
            imin = i;
        }
        if (poly[i].y > ymax) {
            ymax = poly[i].y;
            imax = i;
        }
    }
    double dmax = ymax - ymin;
    double d = dmax;
    int is = imin, js = imax;
    int i = is, j = js;
    do {
        if (cross(poly[(i+1)%n] - poly[i], poly[(j+1)%n] - poly[j]) < 0) {
            i = (i+1)%n;
        } else {
            j = (j+1)%n;
        }
        d = (poly[i] - poly[j]).len();
        if (d > dmax) {
            dmax = d;
            imin = i;
            imax = j;
        }
    } while (i != is || j != js);

    return dmax;
}

// true : divided
// poly1 : left side of l.v, poly2:right side of l.v
bool convex_cut(const Polygon& poly, const Line& l, double eps, Polygon& poly1, Polygon& poly2)
{
    AABB mmb(poly);
    double diag = (mmb.max - mmb.min).len();
    Point center = (mmb.min + mmb.max) * 0.5;

    double c = 0;
    Point pc = project(l, center, &c);
    Point p1 = getpoint(l, c - diag);
    Point p2 = getpoint(l, c + diag);
    Segment line(p1, p2);

    Polygon poly_div;
    vector<int> intpos;
    int n = (int)poly.size();
    for (int i = 0; i < n; i++) {
        double t, dist;
        poly_div.push_back(poly[i]);
        Point h = project(line, poly[i], &t, &dist);
        if (dist < eps) {
            intpos.push_back((int)poly_div.size()-1);
        }
        else {
            Segment edge(poly[i], poly[(i+1)%n]);
            Point x;
            bool is_int = intersect(line, edge, eps, &x);
            if (is_int && !x.equal(poly[(i+1)%n])) {
                intpos.push_back((int)poly_div.size());
                poly_div.push_back(x);
            }
        }
    }
    n = (int)poly_div.size();

    if (intpos.size() < 2) {
        if (ccw(p1, p2, center, eps) >= 0) {
            poly1 = poly;
        } else {
            poly2 = poly;
        }
        return false;
    }
    int x0 = intpos[0];
    int x1 = intpos[1];

    poly1.clear();
    for (int i = 0; i < x0; i++) {
        poly1.push_back(poly_div[i]);
    }
    poly1.push_back(poly_div[x0]);
    poly1.push_back(poly_div[x1]);
    for (int i = x1+1; i < n; i++) {
        poly1.push_back(poly_div[i]);
    }

    poly2.clear();
    for (int i = x0; i <= x1; i++) {
        poly2.push_back(poly_div[i]);
    }

    if ( dot(poly_div[x1] - poly_div[x0], l.v) < 0 ) {
        poly2.swap(poly1);
    }

    return true;
}

#define CGL_4_C
//------------------------------------------------------------------
int main()
{
    const double eps = 1e-10;
    cout << fixed << setprecision(10);

    int n;
    cin >> n;
    Polygon poly;
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        Point p(x, y);
        poly.push_back(p);
    }

    int q;
    cin >> q;
    for (int i = 0 ; i < q; i++) {
        double x1, y1, x2, y2;
        cin >> x1 >> y1 >> x2 >> y2;
        Point p1(x1, y1), p2(x2, y2);
        Polygon poly1, poly2;
        bool cut = convex_cut(poly, Line(p1, (p2-p1).norm()), eps, poly1, poly2);
        // assert(cut);
        cout << area(poly1) << endl;
    }

    return 0;
}
#ifdef CGL_4_B
int main()
{
    const double eps = 1e-10;
    
    int n;
    cin >> n;

    Polygon poly;
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        Point p(x, y);
        poly.push_back(p);
    }

    cout << fixed << setprecision(10);

    cout << convex_diameter(poly) << endl;

    return 0;
}
#endif

#ifdef CGL_4_A
int main()
{
    const double eps = 1e-10;
    
    int n;
    cin >> n;

    Polygon poly;
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        Point p(x, y);
        poly.push_back(p);
    }

    Polygon hull = convex_hull(poly, eps);

    int m = (int)hull.size();

    cout << fixed << setprecision(0);

    cout << m << endl;
    for (int i = 0; i < m; i++) {
        cout << hull[i].x << " " << hull[i].y << endl;
    }
    return 0;
}
#endif
#ifdef CGL_3_C
int main()
{
    int n;
    cin >> n;

    cout << fixed << setprecision(1);
    const double eps = 1e-10;

    Polygon poly;
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        Point p(x, y);
        poly.push_back(p);
    }

    int q;
    cin >> q;

    for (int i = 0; i < q; i++) {
        double x, y;
        cin >> x >> y;
        Point p(x, y);
        cout << contain(poly, p, eps) << endl;
    }

    return 0;
}
#endif

#ifdef CGL_3_B
int main()
{
    int n;
    cin >> n;

    cout << fixed << setprecision(1);
    const double eps = 1e-10;

    Polygon poly;
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        Point p(x, y);
        poly.push_back(p);
    }

    for (int i = 0; i < n; i++) {
        int prev = (i + n -1) % n;
        int next = (i + 1) % n;
        if (ccw(poly[prev], poly[i], poly[next], eps) < 0) {
            cout << 0 << endl;
            return 0;
        }
    }

    cout << 1 << endl;

    return 0;
}
#endif

#ifdef CGL_3_A
int main()
{
    int n;
    cin >> n;

    cout << fixed << setprecision(1);
    const double eps = 1e-10;

    Polygon poly;

    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        Point p(x, y);
        poly.push_back(p);
    }

    double s = 0;
    for (int i = 1; i < n-1; i++) {
        s += 0.5 * signedArea(poly[0], poly[i], poly[i+1]);
    }

    cout << s << endl;

    return 0;
}
#endif

#ifdef CGL_2_D
int main()
{
    int q;
    cin >> q;

    cout << fixed << setprecision(10);
    const double eps = 1e-10;

    for (int i = 0; i < q; i++) {
        double x0, y0, x1, y1, x2, y2, x3, y3;
        cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;

        Point p0(x0, y0);
        Point p1(x1, y1);
        Point p2(x2, y2);
        Point p3(x3, y3);

        Segment s1(p0, p1), s2(p2, p3);
        double d;
        closest(s1, s2, eps, &d);
        cout << d << endl;
    }

    return 0;
}
#endif

#ifdef CGL_2_C

int main()
{
    int q;
    cin >> q;

    cout << fixed << setprecision(8);
    const double eps = 1e-10;

    for (int i = 0; i < q; i++) {
        double x0, y0, x1, y1, x2, y2, x3, y3;
        cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;

        Point p0(x0, y0);
        Point p1(x1, y1);
        Point p2(x2, y2);
        Point p3(x3, y3);

        Segment s1(p0, p1), s2(p2, p3);
        Point x;
        if (intersect(s1, s2, eps, &x)) {
            cout << x.x << " " << x.y << endl;
        }
        else {
            cout << "No intersection!" << endl;
            return 1;
        }
    }

    return 0;
}

#endif
#ifdef CGL_2_B

int main()
{
    int q;
    cin >> q;

    cout << fixed << setprecision(8);
    const double eps = 1e-10;

    for (int i = 0; i < q; i++) {
        double x0, y0, x1, y1, x2, y2, x3, y3;
        cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;

        Point p0(x0, y0);
        Point p1(x1, y1);
        Point p2(x2, y2);
        Point p3(x3, y3);

        Segment s1(p0, p1), s2(p2, p3);

        if (intersect(s1, s2, eps)) {
            cout << 1 << endl;
        }
        else {
            cout << 0 << endl;
        }
    }

    return 0;
}

#endif
#ifdef CGL_2_A

int main()
{
    int q;
    cin >> q;

    cout << fixed << setprecision(8);
    const double eps = 1e-10;

    for (int i = 0; i < q; i++) {
        double x0, y0, x1, y1, x2, y2, x3, y3;
        cin >> x0 >> y0 >> x1 >> y1 >> x2 >> y2 >> x3 >> y3;

        Point p0(x0, y0);
        Point p1(x1, y1);
        Point p2(x2, y2);
        Point p3(x3, y3);

        Vector v0 = (p1-p0).norm();
        Vector v1 = (p3-p2).norm();

        double c = cross(v0, v1);
        double d = dot(v0, v1);

        if (fabs(d) < eps) {
            cout << 1 << endl;
        }
        else if (fabs(c) < eps) {
            cout << 2 << endl;
        }
        else {
            cout << 0 << endl;
        }
    }

    return 0;
}

#endif
//------------------------------------------------------------------
#ifdef CGL_1_C

int main()
{
    double x1, y1, x2, y2;
    cin >> x1 >> y1 >> x2 >> y2;

    Point p1(x1, y1);
    Point p2(x2, y2);
    Line down(p1, (p2-p1).norm());

    double s = 0;
    double t = param(down, p2);

    int q;
    cin >> q;

    cout << fixed << setprecision(8);
    const double eps = 1e-10;

    for (int i = 0; i < q; i++) {
        double x, y;
        cin >> x >> y;
        Point p3(x, y);
        double m, dist;
        Point h = project(down, p3, &m, &dist);
        if (dist < eps) {
            if (m < s) {
                cout << "ONLINE_BACK" << endl;
            }
            else if (s <= m && m <= t) {
                cout << "ON_SEGMENT" << endl;
            }
            else {
                cout << "ONLINE_FRONT" << endl;
            }
        }
        else {
            if (ccw(p1, p2, p3, eps) > 0) {
                cout << "COUNTER_CLOCKWISE" << endl;
            } else {
                cout << "CLOCKWISE" << endl;
            }
        }
    }

    return 0;
}

#endif

#ifdef CGL_1_B

int main()
{
    double x1, y1, x2, y2;
    cin >> x1 >> y1 >> x2 >> y2;

    Point p1(x1, y1);
    Point p2(x2, y2);
    Line down(p1, (p2-p1).norm());

    int q;
    cin >> q;

    cout << fixed << setprecision(8);

    for (int i = 0; i < q; i++) {
        double x, y;
        cin >> x >> y;
        Point q(x, y);
        Point h = symmetry(down, q);

        cout << h.x << " " << h.y << endl;
    }

    return 0;
}

#endif

#ifdef CGL_1_A

int main()
{
    double x1, y1, x2, y2;
    cin >> x1 >> y1 >> x2 >> y2;

    Point p1(x1, y1);
    Point p2(x2, y2);
    Line down(p1, (p2-p1).norm());

    int q;
    cin >> q;

    cout << fixed << setprecision(8);

    for (int i = 0; i < q; i++) {
        double x, y;
        cin >> x >> y;
        Point q(x, y);
        Point h = project(down, q);

        cout << h.x << " " << h.y << endl;
    }

    return 0;
}

#endif

//------------------------------------------------------------------
const long iSeed = 12345;
mt19937 mt(iSeed);

static bool st_GetSampleInteger2( const long numgrid[2], long numpnt, vector<vector<long>>& vviPos )
{
    assert( numpnt > 0 );
    vector<vector<long>> vviPos0;
    set<vector<long>> z;
    for(long i=0; i<numpnt*5; i++) {
        if((long)vviPos0.size()==numpnt) {
            break;
        }
        vector<long> v;
        for(long k=0; k<2; k++) {
            assert( numgrid[k] > 0 );
            long tmp = mt() % numgrid[k];
            v.push_back(tmp);
        }
        if(z.insert(v).second) {
            vviPos0.push_back( v ); 
        }
    }
    if((long)vviPos0.size() != numpnt) {
        return false;
    }
    vviPos.swap( vviPos0 );
    return true;
}

#ifdef TEST_GENERATOR
int main()
{
    long numgrid[2] = {20000, 2000};
    long numpnt = 10000;

    vector<vector<long>> vviPos;
    if (!st_GetSampleInteger2(numgrid, numpnt, vviPos)) {
        return 1;
    }

    ofstream out("E:\\Develop2\\AoJ\\Courses\\convex_hull.txt");
    if (!out) {
        return 2;
    }

    out << numpnt << endl;
    for (int i = 0; i < numpnt; i++) {
        out << vviPos[i][0] - 10000 << " " << vviPos[i][1] - 10000 << endl;
    }

    return 0;
}
#endif