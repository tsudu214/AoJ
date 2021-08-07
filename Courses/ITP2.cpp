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

using namespace std;

int main()
{
    int n;
    cin >> n;
    int k;
    cin >> k;

    int max = 1;
    int shift = 0;
    for (int i = 0; i < max; i++ ) {
        bitset<32> b(i);
        if ( b.count() == k ) {
            cout << i << ":";
            for (int e = 0; e < n; e++ ) {
                if ( b.test(e) ) {
                    cout << " " << e;
                }
            }
            cout << endl;
        }
        if (i == max-1 && shift < n) {
            max <<= 1;
            shift++;
        }
    }

    return 0;
}

#if 0

int main()
{
    int n;
    cin >> n;

    bitset<32> m(0);
    int k;
    cin >> k;
    for (int j = 0; j < k; j++) {
        int f;
        cin >> f;
        m.set(f);
    }

    int pow2n = (int)pow(2, n);
    for (int i = 0; i < pow2n; i++ ) {
        bitset<32> b(i);
        if ( (b & ~m).none() ) {
            cout << i << ":";
            for (int e = 0; e < n; e++ ) {
                if ( b.test(e) ) {
                    cout << " " << e;
                }
            }
            cout << endl;
        }
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    bitset<32> m(0);
    int k;
    cin >> k;
    for (int j = 0; j < k; j++) {
        int f;
        cin >> f;
        m.set(f);
    }

    int pow2n = pow(2, n);
    for (int i = 0; i < pow2n; i++ ) {
        bitset<32> b(i);
        if ( (b | ~m).all() ) {
            cout << i << ":";
            for (int e = 0; e < n; e++ ) {
                if ( b.test(e) ) {
                    cout << " " << e;
                }
            }
            cout << endl;
        }
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    int m = 1; // mask
    int pow2n = pow(2, n);
    for (int i = 0; i < pow2n; i++ ) {
        cout << i << ":";
        bitset<32> b(i);
        for (int e = 0; e < n; e++ ) {
            if ( b.test(e) ) {
                cout << " " << e;
            }
        }
        cout << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    bitset<64> b(0);

    vector<bitset<64>> m(n, 0);
    for (int i = 0; i < n; i++ ) {
        int k;
        cin >> k;
        for (int j = 0; j < k; j++) {
            int f;
            cin >> f;
            m[i].set(f);
        }
    }

    int q;
    cin >> q;
    for (int iq = 0; iq < q; iq++ ) {
        int c, i;
        cin >> c >> i;
        switch (c) {
        case 0:
            cout << (b.test(i)? 1 : 0) << endl;
            break;
        case 1:
            b |= m[i];
            break;
        case 2:
            b &= (~m[i]);
            break;
        case 3:
            b ^= m[i];
            break;
        case 4:
            cout << ((b | ~m[i]).all()? 1 : 0) << endl;
            break;
        case 5:
            cout << ((b & m[i]).any()? 1 : 0) << endl;
            break;
        case 6:
            cout << ((b & m[i]).none()? 1 : 0) << endl;
            break;
        case 7:
            cout << (b & m[i]).count() << endl;
            break;
        case 8:
            cout << (b & m[i]).to_ullong() << endl;
            break;
        }
    }

    return 0;
}

int main()
{
    bitset<64> b(0);

    int q;
    cin >> q;
    for (int iq = 0; iq < q; iq++ ) {
        int c, i;
        cin >> c;
        switch (c) {
        case 0:
            cin >> i;
            cout << (b.test(i)? 1 : 0) << endl;
            break;
        case 1:
            cin >> i;
            b.set(i);
            break;
        case 2:
            cin >> i;
            b.reset(i);
            break;
        case 3:
            cin >> i;
            b.flip(i);
            break;
        case 4:
            cout << (b.all()? 1 : 0) << endl;
            break;
        case 5:
            cout << (b.any()? 1 : 0) << endl;
            break;
        case 6:
            cout << (b.none()? 1 : 0) << endl;
            break;
        case 7:
            cout << b.count() << endl;
            break;
        case 8:
            cout << b.to_ullong() << endl;
            break;
        }
    }

    return 0;
}

int main()
{
    uint32_t Na, Nb;
    cin >> Na >> Nb;

    bitset<32> a(Na), b(Nb);

    cout << (a & b) << endl;

    cout << (a | b) << endl;

    cout << (a ^ b) << endl;

    return 0;
}

int main()
{
    uint32_t x;
    cin >> x;

    bitset<32> b(x);

    cout << b << endl;

    bitset<32> c = b;

    c.flip();
    cout << c << endl;

    c = b;
    c <<= 1;
    cout << c << endl;

    c = b;
    c >>= 1;
    cout << c << endl;

    return 0;
}

int main()
{
    set<long long> S1, S2;

    int n;
    long long x;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> x;
        S1.insert(x);
    }

    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> x;
        S2.insert(x);
    }

    vector<long long> v;
    set_symmetric_difference(S1.begin(), S1.end(), S2.begin(), S2.end(), back_inserter(v));

    for ( auto& e : v ) {
        cout << e << endl;
    }

    return 0;
}

int main()
{
    set<long long> S1, S2;

    int n;
    long long x;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> x;
        S1.insert(x);
    }

    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> x;
        S2.insert(x);
    }

    vector<long long> v;
    set_difference(S1.begin(), S1.end(), S2.begin(), S2.end(), back_inserter(v));

    for ( auto& e : v ) {
        cout << e << endl;
    }

    return 0;
}

int main()
{
    set<long long> S1, S2;

    int n;
    long long x;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> x;
        S1.insert(x);
    }

    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> x;
        S2.insert(x);
    }

    vector<long long> v;
    set_intersection(S1.begin(), S1.end(), S2.begin(), S2.end(), back_inserter(v));

    for ( auto& e : v ) {
        cout << e << endl;
    }

    return 0;
}

int main()
{
    set<long long> S1, S2;

    int n;
    long long x;
    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> x;
        S1.insert(x);
    }

    cin >> n;
    for (int i = 0; i < n; i++) {
        cin >> x;
        S2.insert(x);
    }

    vector<long long> v;
    set_union(S1.begin(), S1.end(), S2.begin(), S2.end(), back_inserter(v));

    for ( auto& e : v ) {
        cout << e << endl;
    }

    return 0;
}

int main()
{
    multimap<string, long long> M;

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int c;
        string key;
        cin >> c >> key;
        if (c == 0) {
            long long x;
            cin >> x;
            M.insert(make_pair(key, x));
        }
        else if (c == 1) {
            auto range = M.equal_range(key);
            for (auto it = range.first; it != range.second; ++it) {
                cout << it->second << endl;
            }
        }
        else if (c == 2) {
            M.erase(key);
        }
        else if (c == 3) {
            string key2;
            cin >> key2;
            auto itBgn = M.lower_bound(key);
            auto itEnd = M.upper_bound(key2);
            for (auto it = itBgn; it != itEnd; ++it) {
                cout << it->first << " " << it->second << endl;
            }
        }
    }

    return 0;
}

int main()
{
    map<string, long long> M;

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int c;
        string key;
        cin >> c >> key;
        if (c == 0) {
            long long x;
            cin >> x;
            M[key] = x;
        }
        else if (c == 1) {
            if (M.count(key)) {
                cout << M[key] << endl;
            } else {
                cout << 0 << endl;
            }
        }
        else if (c == 2) {
            M.erase(key);
        }
        else if (c == 3) {
            string key2;
            cin >> key2;
            auto itBgn = M.lower_bound(key);
            auto itEnd = M.upper_bound(key2);
            for (auto it = itBgn; it != itEnd; ++it) {
                cout << it->first << " " << it->second << endl;
            }
        }
    }

    return 0;
}

int main()
{
    map<string, long long> M;

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int c;
        string key;
        cin >> c >> key;
        if (c == 0) {
            long long x;
            cin >> x;
            M[key] = x;
        }
        else if (c == 1) {
            if (M.count(key)) {
                cout << M[key] << endl;
            } else {
                cout << 0 << endl;
            }
        }
        else if (c == 2) {
            M.erase(key);
        }
    }

    return 0;
}

int main()
{
    multiset<long long> S;

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int c;
        long long x;
        cin >> c >> x;
        if (c == 0) {
            S.insert(x);
            cout << S.size() << endl;
        }
        else if (c == 1) {
            cout << S.count(x) << endl;
        }
        else if (c == 2) {
            S.erase(x);
        }
        else if (c == 3) {
            long long y;
            cin >> y;
            auto itBgn = S.lower_bound(x);
            auto itEnd = S.upper_bound(y);
            for (auto it = itBgn; it != itEnd; ++it) {
                cout << *it << endl;
            }
        }
    }

    return 0;
}

int main()
{
    set<long long> S;

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int c;
        long long x;
        cin >> c >> x;
        if (c == 0) {
            S.insert(x);
            cout << S.size() << endl;
        }
        else if (c == 1) {
            bool found = S.find(x) != S.end();
            cout << (found? 1 : 0) << endl;
        }
        else if (c == 2) {
            S.erase(x);
        }
        else if (c == 3) {
            long long y;
            cin >> y;
            auto itBgn = S.lower_bound(x);
            auto itEnd = S.upper_bound(y);
            for (auto it = itBgn; it != itEnd; ++it) {
                cout << *it << endl;
            }
        }
    }

    return 0;
}

int main()
{
    set<long long> S;

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int c;
        long long x;
        cin >> c >> x;
        if (c == 0) {
            S.insert(x);
            cout << S.size() << endl;
        }
        else if (c == 1) {
            bool found = S.find(x) != S.end();
            cout << (found? 1 : 0) << endl;
        }
        else if (c == 2) {
            S.erase(x);
        }
    }

    return 0;
}

int main()
{
    set<long long> S;

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int c;
        long long x;
        cin >> c >> x;
        if (c == 0) {
            S.insert(x);
            cout << S.size() << endl;
        }
        else if (c == 1) {
            bool found = S.find(x) != S.end();
            cout << (found? 1 : 0) << endl;
        }
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        long long k;
        cin >> k;

        auto it = equal_range(v.begin(), v.end(), k);
        cout << it.first - v.begin() << " " << it.second - v.begin() << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        long long k;
        cin >> k;

        auto it =lower_bound(v.begin(), v.end(), k);
        cout << it - v.begin() << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v1(n);
    for (int i = 0; i < n; i++) {
        cin >> v1[i];
    }
    cin >> n;
    vector<long long> v2(n);
    for (int i = 0; i < n; i++) {
        cin >> v2[i];
    }

    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());

    if (includes(v1.begin(), v1.end(), v2.begin(), v2.end())) {
        cout << 1 << endl;
    } else {
        cout << 0 << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    sort(v.begin(), v.end());

    int q;
    cin >> q;
    for ( int i = 0; i < q; i++ ) {
        int k;
        cin >> k;

        if (binary_search(v.begin(), v.end(), k)) {
            cout << 1 << endl;
        } else {
            cout << 0 << endl;
        }
    }

    return 0;
}

void print_v(const std::vector<int>& v)
{
    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) cout << " ";
        cout << *it;
    }
    cout << endl;
}

int main()
{
    int n;
    cin >> n;
    vector<int> v(n);
    for (int i = 0; i < n; i++) {
        v[i] = i+1;
    }

    do {
        print_v(v);
    } while ( next_permutation(v.begin(), v.end()));

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<int> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    vector<int> prev = v;
    vector<int> next = v;
    bool b1 = prev_permutation(prev.begin(), prev.end());
    bool b2 = next_permutation(next.begin(), next.end());

    if ( b1 ) print_v(prev);
    print_v(v);
    if ( b2 ) print_v(next);

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<tuple<long long, long long, char, long long, string>> a(n);
    for (int i = 0; i < n; i++) {
        long long v, w, d;
        char t;
        string s;
        cin >> v >> w >> t >> d >> s;
        a[i] = make_tuple(v, w, t, d, s);
    }

    sort(a.begin(), a.end());

    for (auto it = a.begin(); it != a.end(); ++it) {
        cout << get<0>(*it) << " " 
             << get<1>(*it) << " " 
             << get<2>(*it) << " "
             << get<3>(*it) << " "
             << get<4>(*it) <<endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<pair<long long, long long>> v(n);
    for (int i = 0; i < n; i++) {
        long long x, y;
        cin >> x >> y;
        v[i] = make_pair(x, y);
    }

    sort(v.begin(), v.end());

    for (auto it = v.begin(); it != v.end(); ++it) {
        cout << it->first << " " << it->second << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    v.erase(unique(v.begin(), v.end()), v.end());

    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) cout << " ";
        cout << *it;
    }
    cout << endl;

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int b, e, t;
        cin >> b >> e >> t;
        swap_ranges(v.begin() + b, v.begin() + e, v.begin() + t);
    }

    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) cout << " ";
        cout << *it;
    }
    cout << endl;

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int b, m, e;
        cin >> b >> m >> e;
        rotate(v.begin() + b, v.begin() + m, v.begin() + e);
    }

    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) cout << " ";
        cout << *it;
    }
    cout << endl;

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int b, e;
        cin >> b >> e;
        reverse(v.begin() + b, v.begin() + e);
    }

    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) cout << " ";
        cout << *it;
    }
    cout << endl;

    return 0;
}

int main()
{
    int n;
    cin >> n;
    vector<long long> v1(n);
    for (int i = 0; i < n; i++) {
        cin >> v1[i];
    }
    int m;
    cin >> m;
    vector<long long> v2(m);
    for (int i = 0; i < m; i++) {
        cin >> v2[i];
    }

    if (v1 < v2) {
        cout << 1 << endl;
    }
    else {
        cout << 0 << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }
    
    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int b, e;
        long long k;
        cin >> b >> e >> k;
        cout << count(v.begin() + b, v.begin() + e, k) << endl;
    }

    return 0;
}

int main()
{
    int n;
    cin >> n;

    vector<long long> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }
    
    int q;
    cin >> q;
    for (int i = 0; i < q; i++) {
        int c, b, e;
        cin >> c >> b >> e;
        if (c == 0) { // min in [b,e)
            cout << *min_element(v.begin() + b, v.begin() + e) << endl;
        } 
        else if (c == 1) { // max in [b, e)
            cout << *max_element(v.begin() + b, v.begin() + e) << endl;
        }

    }

    return 0;
}

int main()
{
    long long a, b, c;
    cin >> a >> b >> c;

    cout << min({a, b, c}) << " " << max({a, b, c}) << endl;

    return 0;
}

int main()
{
    int n, q;
    cin >> n >> q;

    vector<list<long long>> li(n);
    for (int i = 0; i < q; i++) {
        int c, t, s;
        long long a;
        cin >> c >> t;
        if (c == 0) {
            cin >> a;
            li[t].push_back(a);
        }
        else if (c == 1){
            auto it = li[t].begin();
            for (; it != li[t].end(); ++it) {
                if (it!= li[t].begin()) cout << " ";
                cout << *it;
            }
            cout << endl;
        }
        else if (c == 2) {
            cin >> s;
            li[s].splice(li[s].end(), li[t]);
        }
    }

    return 0;
}

int main()
{
    int n, q;
    cin >> n >> q;

    vector<priority_queue<long long>> Q(n);
    for (int i = 0; i < q; i++) {
        int c, t;
        long long a;
        cin >> c >> t;
        if (c == 0) {
            cin >> a;
            Q[t].push(a);
        }
        else if (c == 1){
            if (!Q[t].empty()) {
                long long b = Q[t].top();
                cout << b << endl;
            }
        }
        else if (c == 2) {
            if (!Q[t].empty()) {
                Q[t].pop();
            }
        }
    }

    return 0;
}

int main()
{
    int n, q;
    cin >> n >> q;

    vector<queue<long long>> Q(n);
    for (int i = 0; i < q; i++) {
        int c, t;
        long long a;
        cin >> c >> t;
        if (c == 0) {
            cin >> a;
            Q[t].push(a);
        }
        else if (c == 1){
            if (!Q[t].empty()) {
                long long b = Q[t].front();
                cout << b << endl;
            }
        }
        else if (c == 2) {
            if (!Q[t].empty()) {
                Q[t].pop();
            }
        }
    }

    return 0;
}

int main()
{
    int n, q;
    cin >> n >> q;

    vector<stack<long long>> s(n);
    for (int i = 0; i < q; i++) {
        int c, t;
        long long a;
        cin >> c >> t;
        if (c == 0) {
            cin >> a;
            s[t].push(a);
        }
        else if (c == 1){
            if (!s[t].empty()) {
                long long b = s[t].top();
                cout << b << endl;
            }
        }
        else if (c == 2) {
            if (!s[t].empty()) {
                s[t].pop();
            }
        }
    }

    return 0;
}

int main()
{
    int n, q;
    cin >> n >> q;

    vector<vector<long long>> vv(n);
    for (int i = 0; i < q; i++) {
        int c, t;
        long long a;
        cin >> c >> t;
        if (c == 0) {
            cin >> a;
            vv[t].push_back(a);
        }
        else if (c == 1){
            for (size_t j = 0; j < vv[t].size(); j++) {
                cout << vv[t][j];
                if (j < vv[t].size()-1) cout << " ";
            }
            cout << endl;
        }
        else if (c == 2) {
            vv[t].clear();
        }
    }

    return 0;
}

int main()
{
    int q;
    cin >> q;

    list<long long> li;
    auto it = li.end();
    for (int i = 0; i < q; i++) {
        int c, d;
        long long a;
        cin >> c;
        if (c == 0) {
            cin >> a;
            it = li.insert(it, a);
        }
        else if (c == 1){
            cin >> d;
            if (d > 0) {
                for (int j = 0; j < d; j++) it++;
            }
            else if (d < 0) {
                for (int j = 0; j < -d; j++) it--;
            }
        }
        else if (c == 2) {
            it = li.erase(it);
        }
    }

    for (const auto& e : li) {
        cout << e << endl;
    }

    return 0;
}

int main()
{
    int q;
    cin >> q;

    deque<long long> dq;
    for (int i = 0; i < q; i++) {
        int c, d;
        long long a;
        int p;
        cin >> c;
        if (c == 0) {
            cin >> d >> a;
            if (d == 0)
                dq.push_front(a);
            else if (d == 1)
                dq.push_back(a);
        }
        else if (c == 1){
            cin >> p;
            cout << dq[p] << endl;
        }
        else if (c == 2) {
            cin >> d;
            if (d == 0)
                dq.pop_front();
            else if (d == 1)
                dq.pop_back();
        }
    }

    return 0;
}

int main()
{
    int q;
    cin >> q;

    vector<long long> v;
    for (int i = 0; i < q; i++) {
        int c;
        long long a;
        int p;
        cin >> c;
        if (c == 0) {
            cin >> a;
            v.push_back(a);
        }
        else if (c == 1){
            cin >> p;
            cout << v[p] << endl;
        }
        else if (c == 2) {
            v.pop_back();
        }
    }

    return 0;
}

#endif

