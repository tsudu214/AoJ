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

using namespace std;

using ll = long long; // 9,223,372,036,854,775,807 = 9*10^18


int main()
{

    return 0;
}

#ifdef Problem8

struct Point {
    Point() = default;
    Point(ll x, ll y) : x(x), y(y) {}
    ll x, y;
};

struct Aabb {
    Aabb() = default;
    Aabb(Point m, Point M) : m(m), M(M) {}
    Point m;
    Point M;
};

int main()
{
    int n;
    cin >> n;

    vector<ll> x, y;
    x.reserve(2*n);
    y.reserve(2*n);

    vector<Aabb>  vRect(n);
    for (int i = 0; i < n; i++) {
        ll x1, x2, y1, y2;
        cin >> x1 >> y1 >> x2 >> y2;
        vRect[i] = Aabb(Point(x1, y1), Point(x2, y2));
        x.push_back(x1);
        x.push_back(x2);
        y.push_back(y1);
        y.push_back(y2);
    }

    sort(x.begin(), x.end());
    x.erase(unique(x.begin(), x.end()), x.end());
    sort(y.begin(), y.end());
    y.erase(unique(y.begin(), y.end()), y.end());

    for (auto & r : vRect) {
        r.m.x = lower_bound(x.begin(), x.end(), r.m.x) - x.begin();
        r.m.y = lower_bound(y.begin(), y.end(), r.m.y) - y.begin();
        r.M.x = lower_bound(x.begin(), x.end(), r.M.x) - x.begin();
        r.M.y = lower_bound(y.begin(), y.end(), r.M.y) - y.begin();
    }

    int xn = (int)x.size();
    int yn = (int)y.size();
    vector<vector<pair<ll, int>>> xevent(yn-1);  // { x, 0(start) or 1(end)}

    for (auto & r : vRect) {
        for (ll iy = r.m.y; iy < r.M.y; iy++) {
            xevent[iy].push_back(make_pair(r.m.x, 0));
            xevent[iy].push_back(make_pair(r.M.x, 1));
        }
    }

    vector<vector<pair<ll, ll>>> xrange(yn-1);  // {start, end}
    for (int i = 0; i < yn-1; i++) {
        auto & event = xevent[i];
        sort(event.begin(), event.end());
        int mult = 0;
        for (int j = 0; j < event.size(); j++) {
            if (event[j].second == 0) {
                mult++;
                if (mult == 1) {
                    xrange[i].push_back(make_pair(event[j].first, -1));
                }
            }
            else {
                mult--;
                if (mult == 0) {
                    xrange[i].back().second = event[j].first;
                }
            }
        }
    }

    ll area = 0;
    for (int i = 0; i < yn-1; i++) {
        ll height = y[i+1] - y[i];
        for ( auto& range : xrange[i] ) {
            ll start = x[range.first];
            ll end   = x[range.second];
            ll width = end - start;
            area += width * height;
        }
    }

    cout << area << endl;

    return 0;
}

#endif

#ifdef Problem2

struct Point {
    double x;
    double y;
};

struct Circle {
    Point o;
    double r;
};

struct Segment {
    int c0;
    int c1;
};

double sdist(Point p, Point q) 
{
    return (p.x - q.x)*(p.x - q.x) + (p.y - q.y)*(p.y - q.y);
}

int main()
{
    int n, m;

    while (cin >> n) {
        cin >> m;
        if (n == 0 && m == 0) break;

        vector<Circle> vC;
        for (int i = 0; i < n; i++) {
            int x, y, r;
            cin >> x >> y >> r;
            Point o = {(double)x, (double)y};
            Circle c = {o, (double)r};
            vC.push_back(c);
        }

        vector<Segment> vS;
        for (int i = 0; i < n; i++) {
            for (int j = i+1; j < n; j++) {
                if ( sdist(vC[i].o, vC[j].o) < (vC[i].r + vC[j].r)*(vC[i].r + vC[j].r) ) {
                    Segment s = {i, j};
                    vS.push_back(s);
                }
            }
        }

        for (int i = 0; i < m; i++) {
            if (i > 0 ) cout << " ";

            int px, py, qx, qy;
            cin >> px >> py >> qx >> qy;
            Point p = {(double)px, (double)py};
            Point q = {(double)qx, (double)qy};

            string sp;
            string sq;
            for (int j = 0; j < n; j++) {
                if (sdist(vC[j].o, p) < vC[j].r*vC[j].r) {
                    sp += "1";
                } else {
                    sp += "0";
                }
            }
            for (int j = 0; j < n; j++) {
                if (sdist(vC[j].o, q) < vC[j].r*vC[j].r) {
                    sq += "1";
                } else {
                    sq += "0";
                }
            }

            if (sp == sq) {
                cout << "YES";
            }
            else {
                cout << "NO";
            }
        }
        cout << endl;
    }

    return 0;
}

#endif

#ifdef Problem7

struct Grid {
    long x;
    long y;
    long z;

    bool operator<(const Grid &r) const {
        return x < r.x || (x == r.x && y < r.y) || (x == r.x && y == r.y && z < r.z);
    }
};

int main()
{
    int n;
    cin >> n;

    set<Grid>  G;
    for (int i = 0; i < n; i++) {
        vector<int> v(3);
        cin >> v[0] >> v[1] >> v[2];
        sort(v.begin(), v.end());

        Grid g = {v[0], v[1], v[2]};
        G.insert(g);
    }

    cout << n - (int)G.size() << endl;

    return 0;
}

#endif

#ifdef Problem3

int w, h;

int id(int i, int j) {
    return i * w + j;
}
void cod(int id, int* i, int* j) 
{
    *i = id / w;
    *j = id - (*i) * w;
}

class DisjointSet
{
public:
    DisjointSet(int size = 1){
        rank.resize(size, 0);
        p.resize(size, 0);
        for (int i = 0; i < size; i++) {
            p[i] = i;
        }
    }

    void unite(int x, int y)
    {
        x = find(x);
        y = find(y);

        // ??????????????????????????
        if (rank[x] > rank[y]) {
            p[y] = x;
        }
        else {
            p[x] = y;
            if (rank[x] == rank[y]) {
                rank[y]++;
            }
        }
    }

    int find(int x) 
    {
        if (p[x] == x) {
            return x;
        }
        else {  // ?o?H???k?iroot???????????????????j
            return p[x] = find(p[x]);
        }
    }

    int count(vector<vector<int>>& C)
    {
        set<int> S;
        for (int id = 0; id < (int)p.size(); id++) {
            int i, j;
            cod(id, &i, &j);
            if (C[i][j] == 1) {
                S.insert(find(id));
            }
        }
        return (int)S.size();
    }

private:
    vector<int> rank;
    vector<int> p;

};

int main()
{

    while (cin >> w) {
        cin >> h;
        if (w == 0 && h == 0) break;
        vector<vector<int>> C;
        DisjointSet DS(w*h);
        for (int i = 0; i < h; i++) {
            vector<int> L(w, -1);
            for (int j = 0; j < w; j++) {
                cin >> L[j];
            }
            C.push_back(L);
        }
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                if (C[i][j] == 0) continue;
                int nx[8] = {-1,  0,  1, -1, 1, -1, 0, 1};
                int ny[8] = {-1, -1, -1,  0, 0,  1, 1, 1};
                for (int k = 0; k < 8; k++) {
                    int x = i + nx[k];
                    int y = j + ny[k];
                    if (x >= 0 && x < h && y >= 0 && y < w) {
                        if (C[i][j] == C[x][y]) {
                            DS.unite(id(i, j), id(x, y));
                        }
                    }
                }
            }
        }
        
        cout << DS.count(C) << endl;
    }

    return 0;
}

#endif

#ifdef Problem4

int main()
{
    int n;
    cin >> n;

    set<int> zero;
    set<int> one;
    for (int i = 0; i < n; i++) {
        zero.insert(i);
    }

    for (int i = 0; i <= 2010; i++) {
        set<int> move0, move1;
        for (int t = 0; t < 3; t++) {
            if (!zero.empty()) {
                auto it = zero.begin();
                move0.insert(*it);
                zero.erase(*it);
            }
            else if (!one.empty()) {
                auto it = one.begin();
                move1.insert(*it);
                one.erase(*it);
            }
        }

        for (auto m : move0) {
            one.insert(m);
        }

        if ( zero.empty() && one.empty() ){
            cout << i+1 << endl;
            return 0;
        }
    }

    return 0;
}

#endif


