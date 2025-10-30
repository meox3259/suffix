#include <iostream>
#include <vector>
#include <climits>
#include <algorithm>

class MaxFenwickTree {
private:
    std::vector<int> tree;
    std::vector<int> arr;  // 原始数组，用于维护实际值
    int n;
    
    // 获取lowbit：x & (-x)，提取最低位的1
    int lowbit(int x) {
        return x & (-x);
    }

public:
    // 构造函数：初始化指定大小的树状数组
    MaxFenwickTree(int size) {
        n = size;
        tree.resize(n + 1, INT_MIN);  // 树状数组从1开始索引
        arr.resize(n + 1, INT_MIN);   // 原始数组也从1开始索引
    }
    
    // 构造函数：使用数组初始化
    MaxFenwickTree(const std::vector<int>& initArr) {
        n = initArr.size();
        tree.resize(n + 1, INT_MIN);
        arr.resize(n + 1, INT_MIN);
        
        // 初始化数组
        for (int i = 0; i < n; i++) {
            updateMax(i + 1, initArr[i]);  // 转换为1-based索引
        }
    }
    
    // 单点最大值更新：如果val更大则更新位置idx
    void updateMax(int idx, int val) {
      if (idx < 1 || idx > n)
        return; // 边界检查

      // 只有当新值更大时才更新
      if (val > arr[idx]) {
        int oldVal = arr[idx];
        arr[idx] = val;

        // 需要重建受影响的树节点
        rebuildAffectedNodes(idx);
      }
    }
    
    // 直接设置单点值
    void setValue(int idx, int val) {
        if (idx < 1 || idx > n) return;
        
        arr[idx] = val;
        rebuildAffectedNodes(idx);
    }
    
    // 查询前缀最大值：从1到idx的最大值
    int queryPrefixMax(int idx) {
        if (idx < 1) return INT_MIN;
        if (idx > n) idx = n;
        
        int maxVal = INT_MIN;
        while (idx > 0) {
            maxVal = std::max(maxVal, tree[idx]);
            idx -= lowbit(idx);
        }
        
        return maxVal;
    }
    
    // 查询区间最大值：从l到r的最大值
    int queryRangeMax(int l, int r) {
        if (l < 1 || r > n || l > r) return INT_MIN;
        
        // 对于最大值查询，不能简单用前缀相减
        // 需要直接遍历区间
        int maxVal = INT_MIN;
        for (int i = l; i <= r; i++) {
            maxVal = std::max(maxVal, arr[i]);
        }
        return maxVal;
    }
    
    // 获取单点值
    int getValue(int idx) {
        if (idx < 1 || idx > n) return INT_MIN;
        return arr[idx];
    }
    
    // 获取全局最大值
    int getGlobalMax() {
        return queryPrefixMax(n);
    }
    
    // 打印原始数组
    void printArray() {
        std::cout << "Array: [";
        for (int i = 1; i <= n; i++) {
            std::cout << arr[i];
            if (i < n) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    // 打印树状数组结构（调试用）
    void printTree() {
        std::cout << "Fenwick Tree: [";
        for (int i = 1; i <= n; i++) {
            std::cout << tree[i];
            if (i < n) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    // 打印所有前缀最大值
    void printPrefixMax() {
        std::cout << "Prefix Max: [";
        for (int i = 1; i <= n; i++) {
            std::cout << queryPrefixMax(i);
            if (i < n) std::cout << ", ";
        }
        std::cout << "]\n";
    }

private:
    // 重建受影响的树节点
    void rebuildAffectedNodes(int idx) {
        // 找到所有需要重建的节点
        std::vector<int> nodesToRebuild;
        
        // 向上找到所有包含idx的节点
        int temp = idx;
        while (temp <= n) {
            nodesToRebuild.push_back(temp);
            temp += lowbit(temp);
        }
        
        // 重建这些节点
        for (int node : nodesToRebuild) {
            rebuildNode(node);
        }
    }
    
    // 重建单个节点的值
    void rebuildNode(int idx) {
        // 计算该节点负责的区间
        int range = lowbit(idx);
        int start = idx - range + 1;
        int end = idx;
        
        // 计算区间内的最大值
        int maxVal = INT_MIN;
        for (int i = start; i <= end; i++) {
            maxVal = std::max(maxVal, arr[i]);
        }
        
        tree[idx] = maxVal;
    }
};

// 优化版本：使用更高效的查询方法
class OptimizedMaxFenwickTree {
private:
    std::vector<int> tree;
    std::vector<int> arr;
    int n;
    
    int lowbit(int x) {
        return x & (-x);
    }
    
    // 计算节点覆盖的区间最大值
    int calculateNodeMax(int idx) {
        int range = lowbit(idx);
        int start = idx - range + 1;
        int maxVal = INT_MIN;
        
        for (int i = start; i <= idx; i++) {
            maxVal = std::max(maxVal, arr[i]);
        }
        
        return maxVal;
    }

public:
    OptimizedMaxFenwickTree(int size) {
        n = size;
        tree.resize(n + 1, INT_MIN);
        arr.resize(n + 1, INT_MIN);
    }
    
    OptimizedMaxFenwickTree(const std::vector<int>& initArr) {
        n = initArr.size();
        tree.resize(n + 1, INT_MIN);
        arr.resize(n + 1, INT_MIN);
        
        // 直接设置原始数组
        for (int i = 0; i < n; i++) {
            arr[i + 1] = initArr[i];
        }
        
        // 构建树状数组
        build();
    }
    
    // 构建整个树状数组
    void build() {
        for (int i = 1; i <= n; i++) {
            tree[i] = calculateNodeMax(i);
        }
    }
    
    void updateMax(int idx, int val) {
        if (idx < 1 || idx > n) return;
        
        if (val > arr[idx]) {
            arr[idx] = val;
            
            // 更新所有受影响的节点
            int temp = idx;
            while (temp <= n) {
                tree[temp] = std::max(tree[temp], val);
                temp += lowbit(temp);
            }
        }
    }
    
    void setValue(int idx, int val) {
        if (idx < 1 || idx > n) return;
        
        arr[idx] = val;
        
        // 重建所有受影响的节点
        int temp = idx;
        while (temp <= n) {
            tree[temp] = calculateNodeMax(temp);
            temp += lowbit(temp);
        }
    }
    
    int queryPrefixMax(int idx) {
        if (idx < 1) return INT_MIN;
        if (idx > n) idx = n;
        
        int maxVal = INT_MIN;
        while (idx > 0) {
            maxVal = std::max(maxVal, tree[idx]);
            idx -= lowbit(idx);
        }
        
        return maxVal;
    }
    
    int getValue(int idx) {
        if (idx < 1 || idx > n) return INT_MIN;
        return arr[idx];
    }
    
    void printArray() {
        std::cout << "Array: [";
        for (int i = 1; i <= n; i++) {
            std::cout << arr[i];
            if (i < n) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    
    void printPrefixMax() {
        std::cout << "Prefix Max: [";
        for (int i = 1; i <= n; i++) {
            std::cout << queryPrefixMax(i);
            if (i < n) std::cout << ", ";
        }
        std::cout << "]\n";
    }
};
