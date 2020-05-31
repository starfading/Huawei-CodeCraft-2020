#include <bits/stdc++.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
using namespace std;

//#define TEST

#define Inline __inline__ __attribute__((always_inline))

typedef uint8_t ui8;
typedef uint16_t ui16;
typedef uint32_t ui32;
typedef uint64_t ui64;

/**
 * 定义：
 * 线程数目
 * 映射后节点数目
 * 反哈希数组
 * 多线程读取文件的必要信息
 */
const ui8 THREAD_NUM = 12;
ui32 g_node_cnt;
char* g_buf;
struct Bi{
    ui32 tid;
    ui32 l;
    ui32 r;
};
vector<ui32> g_unhash_id;
vector<ui32> g_node_data[THREAD_NUM];
vector<ui32> g_graph_data[THREAD_NUM];


/**
 * 整型转换为字符型函数
 * @param s
 * @return
 */
Inline ui32 MyAtoi(char* s){
    ui32 ans = 0;
    while(47 < *s && *s < 58){
        ans *= 10;
        ans += (*s - 48);
        s++;
    }
    return ans;
}

/**
 * 多线程读取数据的工作函数
 * @param bi
 */
void ReadJob(Bi& bi){
    char* bg = g_buf + bi.l;
    const char* ed = g_buf + bi.r;
    char tmp[16];
    g_node_data[bi.tid].reserve(1000000);
    g_graph_data[bi.tid].reserve(500000);
    ui32 line[3];
    while(bg < ed){
        for(ui32 & ii : line){
            ui8 index = 0;
            while(*(bg+index) != ',' && *(bg+index) != '\n') index++;
            memset(tmp, ' ', 16);
            memcpy(tmp, bg, 16);
            ii = MyAtoi(tmp);
            bg += (index + 1);
        }
        // 去掉金额为0转账记录
        if(line[2] != 0){
            g_node_data[bi.tid].emplace_back(line[0]);
            g_node_data[bi.tid].emplace_back(line[1]);
            g_graph_data[bi.tid].emplace_back(line[2]);
        }
    }
}

/**
 * 多线程读取数据的任务分配函数
 * @param file_name
 */
void ReadData(const string& file_name){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    int fd = open(file_name.c_str(), O_RDONLY);
    ui32 length = lseek(fd, 0, SEEK_END);
    g_buf = (char*)mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
    thread threads[THREAD_NUM];
    Bi bi[THREAD_NUM];
    ui32 pre_p = 0;
    for(ui8 i = 0; i < THREAD_NUM; ++i){
        ui32 p = (i + 1) * length / THREAD_NUM;
        while (p > 0 && g_buf[p - 1] != '\n') --p;
        bi[i].l = pre_p;
        bi[i].r = p;
        pre_p = p;
        bi[i].tid = i;
        threads[i] = thread(ReadJob, ref(bi[i]));
    }
    for(auto & thread : threads){
        thread.join();
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "[read test data time]: " << dr_ms  << " ms\n";
#endif
}


/**
 * 定义：
 * 最大点数目，最大边数目
 * 点的出度，出度开始位置
 * 点的入度，入度开始位置
 * 前向星数据结构
 */
const ui32 MAX_NODE = 2500000;
const ui32 MAX_EDGE = 2500000;
ui32 g_max_edge;
ui64 g_sum_edge;
ui16 g_ind[MAX_NODE];
ui16 g_outd[MAX_NODE];
ui32 g_ind_bg[MAX_NODE];
ui32 g_head[MAX_NODE];
struct Edge {
    ui32 to;
    ui32 money;
};
Edge g_graph[MAX_EDGE];


/**
 * 建图函数：
 * 合并线程读取的数据，排序、去重
 * 获得映射后节点数目，反哈希数组
 * 哈希原始节点
 * 使用前向星建图
 */
void BuildGraph(){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    // 数据处理
    vector<ui32> node_d;
    for(auto & it : g_node_data){
        node_d.insert(node_d.end(), it.begin(), it.end());
    }
    vector<ui32> tmp = node_d;
    sort(tmp.begin(), tmp.end());
    tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
    g_node_cnt = tmp.size();
    g_unhash_id = tmp;
    unordered_map<ui32, ui32> hash_id;
    for(ui32 i = 0; i < g_node_cnt; ++i)
        hash_id[tmp[i]] = i;
    // 开始建图
    ui32* curlen = new ui32[g_node_cnt + 1]();
    for(auto & it : g_node_data){
        for(ui32 j = 0; j < it.size(); j+=2){
            g_outd[hash_id[it[j]]]++;
        }
    }
    g_head[0] = 0;
    for(ui32 i = 1; i <= g_node_cnt; ++i){
        g_head[i] = g_head[i - 1] + g_outd[i - 1];
    }
    ui32 from, to, money;
    for(ui8 i = 0; i < THREAD_NUM; ++i){
        for(ui32 j = 0, k = 0; j < g_node_data[i].size(); j+=2, ++k){
            from = hash_id[g_node_data[i][j]];
            to = hash_id[g_node_data[i][j+1]];
            money = g_graph_data[i][k];
            g_graph[g_head[from] + curlen[from]].to = to;
            g_graph[g_head[from] + curlen[from]++].money = money;
            g_ind[to]++;
            g_max_edge = g_max_edge < money ? money : g_max_edge;
            g_sum_edge += money;
        }
    }
    g_ind_bg[0] = 0;
    for(ui32 i = 1; i < g_node_cnt; ++i){
        g_ind_bg[i] = g_ind_bg[i - 1] + g_ind[i-1];
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "[build graph time]: " << dr_ms  << " ms\n";
#endif
}


/**
 * 定义：
 * 不使用Dijk算法的节点标记
 * 前驱链中出度为1的链路数目总和
 * 前驱数组：当“前驱点”的“前驱链中出度为1的链路数目总和”大于0时，记录
 */
bool g_jump_over[MAX_NODE];
ui32 g_pre_single_outd_num[MAX_NODE];
vector<ui32> g_pre_node[MAX_NODE];


/**
 * 拓扑排序函数：
 * 当拓扑到的点出度大于1时，直接拓扑排序即可
 * 当拓扑到的点出度为1时，标记其不使用Dijk算法，并酌情记录其他信息，继续拓扑排序
 */
void TopoSort(){
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    queue<ui32> q;
    for(int i = 0; i < g_node_cnt; ++i){
        if(0 == g_ind[i]) q.push(i);
    }
    while(!q.empty()){
        int v = q.front();
        q.pop();
        if(g_outd[v] > 1){
            const ui32 &l = g_head[v], &r = g_head[v+1];
            Edge* e = &g_graph[l];
            for(ui32 k = l; k < r; ++k, ++e) if(0 == --g_ind[e->to]) q.push(e->to);
        }
        else if(g_outd[v] == 1){
            g_jump_over[v] = true;
            const ui32 &adj_node = g_graph[g_head[v]].to;
            g_pre_single_outd_num[adj_node] += (g_pre_single_outd_num[v] + 1);
            if(g_pre_single_outd_num[v] > 0) g_pre_node[adj_node].emplace_back(v);
            if(0 == --g_ind[adj_node]) q.push(adj_node);
        }
    }
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double,std::milli>(end-start).count();
    cout << "[topo sort time]: " << dr_ms  << " ms\n";
#endif
}


/**
 * 定义：
 * [速度测试]总处理节点数目--原子量
 * 线程抢占的节点--原子量
 * 线程的关键中心性数组
 */
#ifdef TEST
atomic<ui32> cnt(0);
#endif
atomic<ui32> g_th_node(0);
double* g_bc[THREAD_NUM];


/**
 * 更新没有使用Dijk算法的点的关键中心性
 * 计算公式：对于当前点的前驱数组中记录的点u：g_bc[u] += g_pre_single_outd_num[u] * new_delta;
 * @param s
 * @param delta
 * @param tid
 */
void UpdataJumpOverNode(const ui32 &s, const double &new_delta, const ui8 &tid){
    for(ui32 & u : g_pre_node[s]){
        g_bc[tid][u] += g_pre_single_outd_num[u] * new_delta;
        UpdataJumpOverNode(u, new_delta + 1, tid);
    }
}


/**
 * class Use16ForDst：
 * 使用uint16记录Dijk的遍历距离
 */
class Use16ForDst{
public:
    struct Heap{
        ui32 num;
        ui16 dist;
        bool operator < (const Heap &t) const {
            return dist > t.dist;
        }
    };

public:
    static void DijkWithHeap(const ui32& s, const ui8& tid, ui32* pred, ui16* len_of_pred, ui16* dst, ui16* sigma, double* delta, ui32* lst, bool* visit, Heap* heap) {
        int len_of_lst = 0;
        sigma[s] = 1;

        // Dijk算法 + 堆优化
        Heap tmp{};
        tmp.num = s;
        tmp.dist = 0;
        heap[0] = tmp;
        ui32 heap_sz = 1;
        while(heap_sz) {
            pop_heap(heap, heap + heap_sz--);
            if(visit[heap[heap_sz].num]) continue;
            Heap u = heap[heap_sz];
            const ui32 &v = u.num;
            visit[v] = true;
            lst[len_of_lst++] = v;
            const ui32 &l = g_head[v], &r = g_head[v+1];
            Edge* e = &g_graph[l];
            for(ui32 k = l; k < r; ++k, ++e){
                const ui32 &adj_node = e->to;
                const ui16 &new_dist = u.dist + e->money;
                if(new_dist > dst[adj_node]) continue;
                if(!visit[adj_node] && new_dist < dst[adj_node]){
                    dst[adj_node] = new_dist;
                    tmp.num = adj_node;
                    tmp.dist = new_dist;
                    heap[heap_sz++] = tmp;
                    push_heap(heap, heap + heap_sz);
                    len_of_pred[adj_node] = 0;
                    pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
                    sigma[adj_node] = sigma[v];
                }
                else if(new_dist == dst[adj_node]) {
                    pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
                    sigma[adj_node] += sigma[v];
                }
            }
        }

        // 更新本轮Dijk的关键中心性贡献值
        for(int i = len_of_lst - 1; i > 0; --i) {
            const ui32 &v = lst[i];
            for(ui32 j = 0; j < len_of_pred[v]; ++j){
                delta[pred[g_ind_bg[v] + j]] += (1.0 + delta[v]) * sigma[pred[g_ind_bg[v] + j]] / sigma[v];
            }
            g_bc[tid][v] += ((g_pre_single_outd_num[s] + 1) * delta[v]);
        }
        g_bc[tid][s] += g_pre_single_outd_num[s] * delta[s];
        if(!g_pre_node[s].empty()) UpdataJumpOverNode(s, delta[s] + 1, tid);

        // 为下一次Dijk初始化
        for(ui32 i = 0; i < len_of_lst; ++i) {
            const ui32 &v = lst[i];
            dst[v] = INT16_MAX;
            sigma[v] = 0;
            delta[v] = 0.0;
            visit[v] = false;
        }
    }

    /**
     * 多线程Dijk的工作函数
     * 使用原子量抢占获得Dijk的起始点
     * [测试]每处理1000个点，打印一下
     * @param tid
     */
    static void GetBc(const ui8 tid) {
        bool* visit = new bool[g_node_cnt]();
        ui16* sigma = new ui16[g_node_cnt]();
        ui16* len_of_pred = new ui16[g_node_cnt];
        ui32* pred = new ui32[MAX_EDGE];
        ui32* lst = new ui32[g_node_cnt];
        ui16* dst = new ui16[g_node_cnt];
        auto* delta = new double[g_node_cnt]();
        Heap* heap = new Heap[MAX_NODE];
        for(ui32 i = 0; i < g_node_cnt; ++i){
            dst[i] = INT16_MAX;
        }
        while(true){
            ui32 node = g_th_node++;
            if(node >= g_node_cnt) break;
            if(g_jump_over[node]) continue;
            DijkWithHeap(node, tid, pred, len_of_pred, dst, sigma, delta, lst, visit, heap);
#ifdef TEST
            ui32 t = ++cnt;
            if(t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
        }
        delete [] visit;
        delete [] sigma;
        delete [] len_of_pred;
        delete [] pred;
        delete [] lst;
        delete [] dst;
        delete [] delta;
        delete [] heap;
    }

    /**
     * 多线程Dijk任务分配函数
     */
    static void AllocTask(){
#ifdef TEST
        auto start = std::chrono::steady_clock::now();
#endif
        thread threads[THREAD_NUM];
        for(ui8 i = 0; i < THREAD_NUM; ++i){
            g_bc[i] = new double[g_node_cnt]();
            threads[i] = thread(Use16ForDst::GetBc, i);
        }
        for(auto & thread : threads){
            thread.join();
        }
#ifdef TEST
        auto end = std::chrono::steady_clock::now();
        double dr_ms = std::chrono::duration<double, std::milli>(end - start).count();
        cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
    }
};

/**
 * class Use32ForDst：
 * 使用uint32记录Dijk的遍历距离
 */
class Use32ForDst{
public:
    struct Heap{
        ui32 num;
        ui32 dist;
        bool operator < (const Heap &t) const {
            return dist > t.dist;
        }
    };

public:
    static void DijkWithHeap(const ui32& s, const ui8& tid, ui32* pred, ui16* len_of_pred, ui32* dst, ui16* sigma, double* delta, ui32* lst, bool* visit, Heap* heap) {
        int len_of_lst = 0;
        sigma[s] = 1;

        // Dijk算法 + 堆优化
        Heap tmp{};
        tmp.num = s;
        tmp.dist = 0;
        heap[0] = tmp;
        ui32 heap_sz = 1;
        while(heap_sz) {
            pop_heap(heap, heap + heap_sz--);
            if(visit[heap[heap_sz].num]) continue;
            Heap u = heap[heap_sz];
            const ui32 &v = u.num;
            visit[v] = true;
            lst[len_of_lst++] = v;
            const ui32 &l = g_head[v], &r = g_head[v+1];
            Edge* e = &g_graph[l];
            for(ui32 k = l; k < r; ++k, ++e){
                const ui32 &adj_node = e->to;
                const ui32 &new_dist = u.dist + e->money;
                if(new_dist > dst[adj_node]) continue;
                if(!visit[adj_node] && new_dist < dst[adj_node]){
                    dst[adj_node] = new_dist;
                    tmp.num = adj_node;
                    tmp.dist = new_dist;
                    heap[heap_sz++] = tmp;
                    push_heap(heap, heap + heap_sz);
                    len_of_pred[adj_node] = 0;
                    pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
                    sigma[adj_node] = sigma[v];
                }
                else if(new_dist == dst[adj_node]) {
                    pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
                    sigma[adj_node] += sigma[v];
                }
            }
        }

        // 更新本轮Dijk的关键中心性贡献值
        for(int i = len_of_lst - 1; i > 0; --i) {
            const ui32 &v = lst[i];
            for(ui32 j = 0; j < len_of_pred[v]; ++j){
                delta[pred[g_ind_bg[v] + j]] += (1.0 + delta[v]) * sigma[pred[g_ind_bg[v] + j]] / sigma[v];
            }
            g_bc[tid][v] += ((g_pre_single_outd_num[s] + 1) * delta[v]);
        }
        g_bc[tid][s] += g_pre_single_outd_num[s] * delta[s];
        if(!g_pre_node[s].empty()) UpdataJumpOverNode(s, delta[s] + 1, tid);

        // 为下一次Dijk初始化
        for(ui32 i = 0; i < len_of_lst; ++i) {
            const ui32 &v = lst[i];
            dst[v] = INT32_MAX;
            sigma[v] = 0;
            delta[v] = 0.0;
            visit[v] = false;
        }
    }

    /**
     * 多线程Dijk的工作函数
     * 使用原子量抢占获得Dijk的起始点
     * [测试]每处理1000个点，打印一下
     * @param tid
     */
    static void GetBc(const ui8 tid) {
        bool* visit = new bool[g_node_cnt]();
        ui16* sigma = new ui16[g_node_cnt]();
        ui16* len_of_pred = new ui16[g_node_cnt];
        ui32* pred = new ui32[MAX_EDGE];
        ui32* lst = new ui32[g_node_cnt];
        ui32* dst = new ui32[g_node_cnt];
        auto* delta = new double[g_node_cnt]();
        Heap* heap = new Heap[MAX_NODE];
        for(ui32 i = 0; i < g_node_cnt; ++i){
            dst[i] = INT32_MAX;
        }
        while(true){
            ui32 node = g_th_node++;
            if(node >= g_node_cnt) break;
            if(g_jump_over[node]) continue;
            DijkWithHeap(node, tid, pred, len_of_pred, dst, sigma, delta, lst, visit, heap);
#ifdef TEST
            ui32 t = ++cnt;
            if(t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
        }
        delete [] visit;
        delete [] sigma;
        delete [] len_of_pred;
        delete [] pred;
        delete [] lst;
        delete [] dst;
        delete [] delta;
        delete [] heap;
    }

    /**
     * 多线程Dijk任务分配函数
     */
    static void AllocTask(){
#ifdef TEST
        auto start = std::chrono::steady_clock::now();
#endif
        thread threads[THREAD_NUM];
        for(ui8 i = 0; i < THREAD_NUM; ++i){
            g_bc[i] = new double[g_node_cnt]();
            threads[i] = thread(Use32ForDst::GetBc, i);
        }
        for(auto & thread : threads){
            thread.join();
        }
#ifdef TEST
        auto end = std::chrono::steady_clock::now();
        double dr_ms = std::chrono::duration<double, std::milli>(end - start).count();
        cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
    }
};

/**
 * class Use64ForDst：
 * 使用uint64记录Dijk的遍历距离
 */
class Use64ForDst{
public:
    struct Heap{
        ui32 num;
        ui64 dist;
        bool operator < (const Heap &t) const {
            return dist > t.dist;
        }
    };

public:
    static void DijkWithHeap(const ui32& s, const ui8& tid, ui32* pred, ui16* len_of_pred, ui64* dst, ui16* sigma, double* delta, ui32* lst, bool* visit, Heap* heap) {
        int len_of_lst = 0;
        sigma[s] = 1;


        Heap tmp{};
        tmp.num = s;
        tmp.dist = 0;
        heap[0] = tmp;
        ui32 heap_sz = 1;
        while(heap_sz) {
            pop_heap(heap, heap + heap_sz--);
            if(visit[heap[heap_sz].num]) continue;
            Heap u = heap[heap_sz];
            const ui32 &v = u.num;
            visit[v] = true;
            lst[len_of_lst++] = v;
            const ui32 &l = g_head[v], &r = g_head[v+1];
            Edge* e = &g_graph[l];
            for(ui32 k = l; k < r; ++k, ++e){
                const ui32 &adj_node = e->to;
                const ui64 &new_dist = u.dist + e->money;
                if(new_dist > dst[adj_node]) continue;
                if(!visit[adj_node] && new_dist < dst[adj_node]){
                    dst[adj_node] = new_dist;
                    tmp.num = adj_node;
                    tmp.dist = new_dist;
                    heap[heap_sz++] = tmp;
                    push_heap(heap, heap + heap_sz);
                    len_of_pred[adj_node] = 0;
                    pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
                    sigma[adj_node] = sigma[v];
                }
                else if(new_dist == dst[adj_node]) {
                    pred[g_ind_bg[adj_node] + len_of_pred[adj_node]++] = v;
                    sigma[adj_node] += sigma[v];
                }
            }
        }

        // 更新本轮Dijk的关键中心性贡献值
        for(int i = len_of_lst - 1; i > 0; --i) {
            const ui32 &v = lst[i];
            for(ui32 j = 0; j < len_of_pred[v]; ++j){
                delta[pred[g_ind_bg[v] + j]] += (1.0 + delta[v]) * sigma[pred[g_ind_bg[v] + j]] / sigma[v];
            }
            g_bc[tid][v] += (g_pre_single_outd_num[s] + 1) * delta[v];
        }
        g_bc[tid][s] += g_pre_single_outd_num[s] * delta[s];
        if(!g_pre_node[s].empty()) UpdataJumpOverNode(s, delta[s] + 1, tid);

        // 为下一次Dijk初始化
        for(ui32 i = 0; i < len_of_lst; ++i) {
            const ui32 &v = lst[i];
            dst[v] = INT64_MAX;
            sigma[v] = 0;
            delta[v] = 0.0;
            visit[v] = false;
        }
    }

    /**
     * 多线程Dijk的工作函数
     * 使用原子量抢占获得Dijk的起始点
     * [测试]每处理1000个点，打印一下
     * @param tid
     */
    static void GetBc(const ui8 tid) {
        bool* visit = new bool[g_node_cnt]();
        ui16* sigma = new ui16[g_node_cnt]();
        ui16* len_of_pred = new ui16[g_node_cnt];
        ui32* pred = new ui32[MAX_EDGE];
        ui32* lst = new ui32[g_node_cnt];
        ui64* dst = new ui64[g_node_cnt];
        auto* delta = new double[g_node_cnt]();
        Heap* heap = new Heap[MAX_NODE];
        for(ui32 i = 0; i < g_node_cnt; ++i){
            dst[i] = INT64_MAX;
        }
        while(true){
            ui32 node = g_th_node++;
            if(node >= g_node_cnt) break;
            if(g_jump_over[node]) continue;
            DijkWithHeap(node, tid, pred, len_of_pred, dst, sigma, delta, lst, visit, heap);
#ifdef TEST
            ui32 t = ++cnt;
            if(t % 1000 == 0 || t == g_node_cnt) printf("%d/%d\n", t, g_node_cnt);
#endif
        }
        delete [] visit;
        delete [] sigma;
        delete [] len_of_pred;
        delete [] pred;
        delete [] lst;
        delete [] dst;
        delete [] delta;
        delete [] heap;
    }

    /**
     * 多线程Dijk任务分配函数
     */
    static void AllocTask(){
#ifdef TEST
        auto start = std::chrono::steady_clock::now();
#endif
        thread threads[THREAD_NUM];
        for(ui8 i = 0; i < THREAD_NUM; ++i){
            g_bc[i] = new double[g_node_cnt]();
            threads[i] = thread(Use64ForDst::GetBc, i);
        }
        for(auto & thread : threads){
            thread.join();
        }
#ifdef TEST
        auto end = std::chrono::steady_clock::now();
        double dr_ms = std::chrono::duration<double, std::milli>(end - start).count();
        cout << "[calculate time]: " << dr_ms << " ms\n";
#endif
    }
};


/**
 * 精度控制函数：用于std::sort
 * @param a
 * @param b
 * @return
 */
Inline bool PrecisionCmp(const pair<double,ui32> &a, const pair<double,ui32> &b){
    if(abs(a.first - b.first) > 0.0001) return a.first > b.first;
    else return a.second < b.second;
}

/**
 * 输出函数：汇总多线程的关键中心性，排序，取top100输出
 * @param result_file
 */
void Write(const string &result_file){
    auto* bc_final = new double[g_node_cnt]();
    for(auto & ii : g_bc){
        for(ui32 j = 0; j < g_node_cnt; ++j){
            bc_final[j] += ii[j];
        }
    }
    auto res = new pair<double, ui32>[g_node_cnt];
    for(ui32 i = 0; i < g_node_cnt; ++i){
        res[i] = make_pair(bc_final[i], i);
    }
    sort(res, res + g_node_cnt, PrecisionCmp);
    ui8 index = g_node_cnt < 100 ? g_node_cnt : 100;
    ofstream write(result_file.c_str(), ios::out);
    for(ui8 i = 0; i < index; ++i){
        write << g_unhash_id[res[i].second] << ',' << fixed << setprecision(3) << res[i].first << '\n';
    }
    write.close();
    delete [] bc_final;
    delete [] res;
}


/**
 * 主函数依次执行：
 * 1、读取函数
 * 2、建图函数
 * 3、拓扑排序函数
 * 4、选择Dijk类
 * 5、输出函数
 * @return
 */
int main() {
    const string test_file = "/data/test_data.txt";
    const string result_file = "/projects/student/result.txt";
#ifdef TEST
    auto start = std::chrono::steady_clock::now();
#endif
    ReadData(test_file);
    BuildGraph();
    TopoSort();
    if(g_max_edge < 2000)
        Use16ForDst::AllocTask();
    else if(g_sum_edge < INT32_MAX)
        Use32ForDst::AllocTask();
    else
        Use64ForDst::AllocTask();
   Write(result_file);
#ifdef TEST
    auto end = std::chrono::steady_clock::now();
    double dr_ms = std::chrono::duration<double, std::milli>(end - start).count();
    cout << "[total cost time]: " << dr_ms << " ms\n";
#endif
    return 0;
}
