#define _USE_MATH_DEFINES
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
#include <functional>

using namespace std;

// https://www.slideshare.net/iwiwi/ss-3578491
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

        // í·Ç¢ï˚ÇçÇÇ¢ÇŸÇ§Ç…Ç¬Ç»Ç∞ÇÈ
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
        else {  // åoòHà≥èkÅirootÇ…íºê⁄Ç¬Ç»Ç¨Ç»Ç®Ç∑Åj
            return p[x] = find(p[x]);
        }
    }

private:
    vector<int> rank;
    vector<int> p;

};

// https://qiita.com/drken/items/cce6fc5c579051e64fab
class WeightedDisjointSet
{
public:
    WeightedDisjointSet(int size = 1, int wIni = 0){
        rank.resize(size, 0);
        diff_w.resize(size, wIni);
        p.resize(size, 0);
        for (int i = 0; i < size; i++) {
            p[i] = i;
        }
    }

    int root(int x) 
    {
        if (p[x] == x) {
            return x;
        }
        int r =  root(p[x]);
        diff_w[x] += diff_w[p[x]];
        p[x] = r;
        return r;
    }

    int weight(int x) {
        root(x);
        return diff_w[x];
    }

    // weight(y) - weight(x) = w Ç∆Ç»ÇÈÇÊÇ§Ç… merge Ç∑ÇÈ
    bool relate(int x, int y, int w)
    {
        // root Ç∆ÇÃèdÇ›ç∑ï™ï‚ê≥
        w += weight(x); 
        w -= weight(y);

        // x, y ÇrootÇ÷
        x = root(x);
        y = root(y);
        if (x == y) {  // already related
            return false;
        }
        
        // rank[x] >= rank[y] Ç…Ç»ÇÈÇÊÇ§Ç…ì¸ÇÍë÷Ç¶
        if (rank[x] < rank[y]) {
            swap(x, y);
            w = -w;
        }
        // y Ç x ÇÃâ∫Ç÷
        if (rank[x] == rank[y]) {
            rank[x]++;
        }
        p[y] = x;
        diff_w[y] = w;

        return true;
    }

    int diff(int x, int y) {
        return weight(y) - weight(x);
    }

private:
    vector<int> rank;
    vector<int> p;
    vector<int> diff_w;

};

// https://www.slideshare.net/iwiwi/ss-3578491
// sqrt(N) decomposition
class SqrtBucket
{
public:
    int inf = INT_MAX; 

    SqrtBucket(int n) {
        int b = 1;
        while(true) {
            if (b*b >= n) break;
            b++;
        }
        A = vector<int>(n, inf);
        B = vector<int>(b, inf);
    }

    // i = 0, 1, ..., n-1
    void update(int i, int x) {
        A[i] = x;
        int b = (int)B.size();
        int bi = i / b;
        B[bi] = inf;
        for (int i = bi*b; i < (bi+1)*b && i < (int)A.size(); i++ ) {
            B[bi] = min( B[bi], A[i] );
        }
    }

    // index-range [s, t)
    int find(int s, int t) {
        int b = (int)B.size();
        int si = s / b;
        int ti = t / b;
        int m = inf;
        for (int i = s; i < (si+1)*b && i < t; i++) {
            m = min( m, A[i] );
        }
        for (int i = si+1; i < ti; i++ ) {
            m = min( m, B[i] );
        }
        for (int i = ti*b; i < t && i >= s; i++ ) {
            m = min( m, A[i] );
        }
        return m;
    }

private:
    vector<int> A;
    vector<int> B;
};

// https://www.slideshare.net/iwiwi/ss-3578491
// segment tree
class SegTreeRmQ
{
public:
    int inf = INT_MAX; 

    SegTreeRmQ(int n) : N(1) {
        while (N < n) {
            N *= 2;
        }
        A = vector<int>(2*N - 1, inf);
    }

    // i = 0, 1, ..., n-1
    void update(int i, int x) {
        i += N - 1;
        A[i] = x;
        while (i > 0) {
            i = (i - 1)/2;
            A[i] = min(A[i*2+1], A[i*2+2]);
        }
    }

    // [s, t) Ç…Ç®ÇØÇÈminÇï‘Ç∑ [l,r)ì‡Ç≈çiÇËçûÇ› 
    int query(int s, int t, int k, int l, int r) {
        if (r <= s || t <= l) {
            return inf;
        }
        if (s <= l && r <= t) {
            return A[k];
        }
        int vl = query(s, t, k*2 + 1, l, (l + r)/2);
        int vr = query(s, t, k*2 + 2, (l + r)/2, r);
        return min(vl, vr);
    }

    // index-range [s, t)
    int find(int s, int t) {
        return query(s, t, 0, 0, N);
    }

private:
    int N;
    vector<int> A;
};

// https://www.creativ.xyz/segment-tree-abstraction-979/
// SegTree generic version
template <class T>
class SegTree
{
public:
    SegTree(int n, T _def, function<T(T, T)> _operation,  function<T(T, T)> _updater)
        : def(_def), operation(_operation), updater(_updater)
    {
        N = 1;
        while (N < n) {
            N *= 2;
        }
        A = vector<T>(2*N - 1, def);
    }

    // i = 0, 1, ..., n-1
    void update(int i, T x) {
        i += N - 1;
        A[i] = updater(A[i], x);
        while (i > 0) {
            i = (i - 1)/2;
            A[i] = operation(A[i*2+1], A[i*2+2]);
        }
    }

    // [s, t) Ç…Ç®ÇØÇÈ <operation> Çï‘Ç∑ [l,r)ì‡Ç≈çiÇËçûÇ› 
    T query(int s, int t, int k, int l, int r) {
        if (r <= s || t <= l) {
            return def;
        }
        if (s <= l && r <= t) {
            return A[k];
        }
        int vl = query(s, t, k*2 + 1, l, (l + r)/2);
        int vr = query(s, t, k*2 + 2, (l + r)/2, r);
        return operation(vl, vr);
    }

    // index-range [s, t)
    int query(int s, int t) {
        return query(s, t, 0, 0, N);
    }

private:
    int N;
    T   def;
    vector<T> A;
    function<T(T, T)> operation;
    function<T(T, T)> updater;
};

//------------------------------------------------------------------
// 2D geometry class
struct Point 
{
    Point() = default;
    Point(double x0, double y0) : x(x0), y(y0) {}

    double x;
    double y;
};

struct PointId : public Point
{
    PointId() = default;
    PointId(int _id, double x0, double y0) : id(_id), Point(x0, y0) {}
    bool operator < (const PointId&p) const { return id < p.id; }

    int id;
};

bool compare_x(const Point &a, const Point &b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

bool compare_y(const Point &a, const Point &b) {
    return a.y < b.y || (a.y == b.y && a.x < b.x);
}

typedef Point Vector;

using Polygon = vector<Point>;

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
    bool contain(Point p) const {
        return (min.x <= p.x && p.x <= max.x && min.y <= p.y && p.y <= max.y);
    }
    Point min;
    Point max;
};

static const int NIL = -1;

struct Node {
    Node() : location(NIL), right(NIL), left(NIL) {}
    int location;
    int right;
    int left;
};

class KDTree {
public:
    KDTree() {}
    KDTree(const vector<Point>& _points) {
        int n = (int)_points.size();
        points.reserve(n);
        for (int i = 0; i < n; i++) {
            points.push_back(PointId(i, _points[i].x, _points[i].y));
        }
        T.reserve(n);
        int root = Make(0, n, 0);  // == 0
    }

    // make tree in range [l, r)
    int Make(int l, int r, int depth)
    {
        if (l >= r) return NIL;
        if (depth % 2 == 0) {
            sort(points.begin() + l, points.begin() + r, compare_x);
        } else {
            sort(points.begin() + l, points.begin() + r, compare_y);
        }
        int mid = (l + r) / 2;
        int t = (int)T.size();
        T.push_back(Node());
        Node& node = T.back();
        node.location = mid;
        node.left  = Make( l, mid, depth + 1 );
        node.right = Make( mid + 1, r, depth + 1 );
        return t;
    }

    void Query(int v, const AABB& aabb, int depth, vector<PointId>& ans)
    {
        PointId Pv = points[T[v].location];

        if (aabb.contain(Pv)) {
            ans.push_back(Pv);
        }
        if (depth % 2 == 0) {
            if (T[v].left != NIL && aabb.min.x <= Pv.x) {
                Query(T[v].left, aabb, depth + 1, ans);
            }
            if (T[v].right != NIL && Pv.x <= aabb.max.x) {
                Query(T[v].right, aabb, depth + 1, ans);
            }
        }
        else {
            if (T[v].left != NIL && aabb.min.y <= Pv.y) {
                Query(T[v].left, aabb, depth + 1, ans);
            }
            if (T[v].right != NIL && Pv.y <= aabb.max.y) {
                Query(T[v].right, aabb, depth + 1, ans);
            }
        }
    }

private:
    vector<PointId> points;
    vector<Node> T;
};

int main()
{
    int n, q;

    cin >> n;
    vector<Point> P(n);
    for (int i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y;
        P[i].x = x; P[i].y = y;
    }

    KDTree T(P);

    cin >> q;
    for (int i = 0; i < q; i++) {
        double xmin, xmax, ymin, ymax;
        cin >> xmin >> xmax >> ymin >> ymax;
        AABB aabb(Point(xmin, ymin), Point(xmax, ymax));

        vector<PointId> ans;
        T.Query(0, aabb, 0, ans);
        sort(ans.begin(), ans.end());
        for (auto& q : ans) {
            cout << q.id << endl;
        }
        cout << endl;
    }

    return 0;
}


#ifdef DSL_2_B_Segment_Tree

int main()
{
    int n, q;
    cin >> n >> q;

    SegTree<int> st(n, 0,
        [](int a, int b){ return a + b; }, 
        [](int a, int b){ return a + b; }
        );

    for (int i = 0; i < q; i++) {
        int com, x, y;
        cin >> com >> x >> y;
        if (com == 0) {
            st.update(x-1, y);
        }
        if (com == 1) {
            cout << st.query(x-1, y) << endl;
        }
    }

    return 0;
}

#endif

#ifdef DSL_2_A_Segment_Tree

int main()
{
    int n, q;
    cin >> n >> q;

    SegTree<int> st(n, INT_MAX, 
        [](int a, int b){ return min(a, b);},
        [](int a, int b){ return b;}
    );

    for (int i = 0; i < q; i++) {
        int com, x, y;
        cin >> com >> x >> y;
        if (com == 0) {
            st.update(x, y);
        }
        if (com == 1) {
            cout << st.query(x, y+1) << endl;
        }
    }

    return 0;
}

#endif

#ifdef DSL_2_A_Segment_Tree_RmQ

int main()
{
    int n, q;
    cin >> n >> q;

    SegTreeRmQ st(n);

    for (int i = 0; i < q; i++) {
        int com, x, y;
        cin >> com >> x >> y;
        if (com == 0) {
            st.update(x, y);
        }
        if (com == 1) {
            cout << st.find(x, y+1) << endl;
        }
    }

    return 0;
}

#endif

#ifdef DSL_2_A_RMQ_Sqrt_Bucket

int main()
{
    int n, q;
    cin >> n >> q;

    SqrtBucket sb(n);

    for (int i = 0; i < q; i++) {
        int com, x, y;
        cin >> com >> x >> y;
        if (com == 0) {
            sb.update(x, y);
        }
        if (com == 1) {
            cout << sb.find(x, y+1) << endl;
        }
    }

    return 0;
}

#endif

#ifdef DSL_1_B_Weighted_Union_Find

int main()
{
    int n, q;
    cin >> n >> q;

    WeightedDisjointSet ds(n);

    for (int i = 0; i < q; i++) {
        int com, x, y, z;
        cin >> com;
        if (com == 0) {
            cin >> x >> y >> z;
            ds.relate(x, y, z);
        }
        if (com == 1) {
            cin >> x >> y;
            if (ds.root(x) != ds.root(y)) {
                cout << "?" << endl;
            }
            else {
                cout << ds.diff(x, y) << endl;
            }
        }
    }

    return 0;
}
#endif


#ifdef DSL_1_A_Union_Find
int main()
{
    int n, q;
    cin >> n >> q;

    DisjointSet ds(n);

    for (int i = 0; i < q; i++) {
        int com, x, y;
        cin >> com >> x >> y;
        if (com == 0) {
            ds.unite(x, y);
        }
        else if (com == 1) {
            if (ds.find(x) == ds.find(y)) {
                cout << 1 << endl;
            } else {
                cout << 0 << endl;
            }
        }
    }

    return 0;
}
#endif

