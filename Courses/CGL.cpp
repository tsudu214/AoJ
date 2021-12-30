#include <iostream>
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

using namespace std;

using ll = long long;
const ll INF = LLONG_MAX;

//------------------------------------------------------------------
// common class
struct Point 
{
    Point(double x0 = 0, double y0 = 0) : x(x0), y(y0) {}

    Point operator+(Point q) const { return Point(x + q.x, y + q.y); }
    Point operator-(Point q) const { return Point(x - q.x, y - q.y); }
    Point operator*(double a) const { return Point(a * x, a * y); }

    Point& operator+=(Point& q) { x += q.x; y += q.y; return *this; }
    Point& operator-=(Point& q) { x -= q.x; y -= q.y; return *this; }
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
        return Point(a*x, a*y);
    }

    double x;
    double y;
};

bool compare_x(const Point &a, const Point &b) {
    return a.x == b.x? a.y < b.y : a.x < b.x;
}

bool compare_y(const Point &a, const Point &b) {
    return a.y == b.y? a.x < b.x : a.y < b.y;
}

typedef Point Vector;

using Polygon = vector<Point>;

double clamp(double x, double min, double max) 
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

double dot(Vector p, Vector q) 
{
    return p.x * q.x + p.y * q.y;
}

// outer-product of 2D vectors is scalar
double cross(Vector p, Vector q) {
    return p.x * q.y - p.y * q.x;
} 

struct Line 
{
    Line(Point p0, Vector v0) : p(p0), v(v0) {}
    Point p;
    Vector v;  // should be normalized
};

struct Segment 
{
    Segment(Point a0, Point b0): a(a0), b(b0) {}
    Point a;
    Point b;
    Line line(double* len = nullptr) const {
        if(len) *len = (b - a).len();
        return Line(a, (b - a).norm());
    }
};

double param(Line l, Point p)
{
    return dot(p - l.p, l.v);
}

Point getpoint(Line l, double t) 
{
    return l.p + l.v * t;
}

Point project(Line l, Point p, double* t = nullptr, double* dist = nullptr) 
{
    double th = param(l, p);
    Point h = getpoint(l, th);
    if (dist) *dist = (p - h).len();
    if (t) *t = th;
    return h;
}

Point symmetry(Line l, Point p)
{
    return p + (project(l, p) - p) * 2;
}

double signedArea(Point p0, Point p1, Point p2)
{
    return cross((p1 - p0), (p2 - p0));
}

// 1:ccw, -1:cw, 0:on-line
int ccw(Point p0, Point p1, Point p2, double eps, double* area = nullptr)
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

bool intersect(Segment s1, Segment s2, double eps, Point* x = nullptr)
{
    Point p1 = s1.a;
    Point p2 = s1.b;
    Point p3 = s2.a;
    Point p4 = s2.b;

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

Point project(Segment s, Point p, double* t = nullptr, double* dist = nullptr) 
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

bool intersect(Line l1, Line l2, double eps, double* sx = nullptr, double* tx = nullptr, Point* x = nullptr)
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

void closest(Segment seg1, Segment seg2, double eps, 
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

void aabb(const Polygon& poly, Point* m, Point* M)
{
    double mx = HUGE_VAL, my = HUGE_VAL;
    double Mx = -HUGE_VAL, My = -HUGE_VAL;

    for (const auto& v : poly) {
        mx = min(mx, v.x);
        my = min(my, v.y);
        Mx = max(Mx, v.x);
        My = max(My, v.y);
    }
    *m = Point(mx, my);
    *M = Point(My, My);
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
int contain(const Polygon& poly, Point p, double eps)
{
    int n = (int)poly.size();

    Point m, M;
    aabb(poly, &m, &M);
    double diag = (M - m).len();

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

Polygon convex_hull(const Polygon& poly, double eps)
{
    vector<Point> s = poly;
    sort(s.begin(), s.end(), compare_x);

    int n = (int)s.size();
    if (n < 3) return s;

    vector<Point> u, l;
    u.reserve(n);
    l.reserve(n);
    u.push_back(s[0]);
    u.push_back(s[1]);
    for (int i = 2; i < n; i++) {
        for (int j = (int)u.size()-1; j >= 1; j--) {
            if (ccw(u[j-1], u[j], s[i], eps) == 1) {
                u.pop_back();
            }
        }
        u.push_back(s[i]);
    }

    l.push_back(s[n-1]);
    l.push_back(s[n-2]);
    for (int i = n-3; i >= 0; i--) {
        for (int j = (int)l.size()-1; j >= 1; j--) {
            if (ccw(l[j-1], l[j], s[i], eps) == 1) {
                l.pop_back();
            }
        }
        l.push_back(s[i]);
    }

    reverse(l.begin(), l.end());
    for (int i = (int)u.size()-2; i >= 1; i--) {
        l.push_back(u[i]);
    }

    return l;
}

#define CGL_4_A
//------------------------------------------------------------------
int main()
{
    int n;
    cin >> n;

    cout << fixed << setprecision(0);
    const double eps = 1e-10;

    Polygon poly;
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        Point p(x, y);
        poly.push_back(p);
    }

    Polygon hull = convex_hull(poly, eps);

    int m = (int)hull.size();
    int idx = 0;
    for (int i = 1; i < m; i++) {
        if (compare_y(hull[i], hull[idx])) {
            idx = i;
        }
    }

    cout << m << endl;
    for (int i = 0; i < m; i++) {
        cout << hull[(idx + i)%m].x << " " << hull[(idx + i)%m].y << endl;
    }

    return 0;
}

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
    Line l(p1, (p2-p1).norm());

    double s = 0;
    double t = param(l, p2);

    int q;
    cin >> q;

    cout << fixed << setprecision(8);
    const double eps = 1e-10;

    for (int i = 0; i < q; i++) {
        double x, y;
        cin >> x >> y;
        Point p3(x, y);
        double m, dist;
        Point h = project(l, p3, &m, &dist);
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
    Line l(p1, (p2-p1).norm());

    int q;
    cin >> q;

    cout << fixed << setprecision(8);

    for (int i = 0; i < q; i++) {
        double x, y;
        cin >> x >> y;
        Point q(x, y);
        Point h = symmetry(l, q);

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
    Line l(p1, (p2-p1).norm());

    int q;
    cin >> q;

    cout << fixed << setprecision(8);

    for (int i = 0; i < q; i++) {
        double x, y;
        cin >> x >> y;
        Point q(x, y);
        Point h = project(l, q);

        cout << h.x << " " << h.y << endl;
    }

    return 0;
}

#endif
