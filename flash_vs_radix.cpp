#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <intrin.h>

// フラッシュソートの実装
void flashSort(std::vector<int>& arr) {
    int n = arr.size();
    if (n <= 1) {
        return; // 要素数が1以下ならソートの必要なし
    }

    // Step 1: 最小値と最大値を見つける
    int minVal = arr[0], maxVal = arr[0];
    for (int i = 1; i < n; ++i) {
        if (arr[i] < minVal) {
            minVal = arr[i];
        } else if (arr[i] > maxVal) {
            maxVal = arr[i];
        }
    }

    // Step 2: バケットの幅を計算
    const int numBuckets = 100000; // バケットの数を指定します（mの値）
    double c = (double)(numBuckets - 1) / (maxVal - minVal);

    // Step 3: 各バケットの要素数をカウントする
    std::vector<int> bucketCounts(numBuckets, 0);
    for (int i = 0; i < n; ++i) {
        int bucketIndex = int(c * (arr[i] - minVal));
        ++bucketCounts[bucketIndex];
    }

    // Step 4: バケットの累積和を計算する
    std::vector<int> prefixSum(numBuckets, 0);
    for (int i = 1; i < numBuckets; ++i) {
        prefixSum[i] = prefixSum[i - 1] + bucketCounts[i - 1];
    }

    // Step 5: バケットに要素を配置する
    std::vector<int> sortedArr(n);
    for (int i = 0; i < n; ++i) {
        int bucketIndex = int(c * (arr[i] - minVal));
        sortedArr[prefixSum[bucketIndex]] = arr[i];
        ++prefixSum[bucketIndex];
    }

    // Step 6: 各バケットを挿入ソートでソートする
    for (int i = 0; i < numBuckets; ++i) {
        int start = (i == 0) ? 0 : prefixSum[i - 1];
        int end = prefixSum[i];
        // 挿入ソート
        for (int j = start + 1; j < end; ++j) {
            int key = sortedArr[j];
            int k = j - 1;
            while (k >= start && sortedArr[k] > key) {
                sortedArr[k + 1] = sortedArr[k];
                k--;
            }
            sortedArr[k + 1] = key;
        }
    }

    // ソート結果を入力配列に戻す
    arr = sortedArr;
}

// 基数ソートの実装
void radixSort(std::vector<int>& arr) {
    // 最大値を見つける
    int maxVal = *std::max_element(arr.begin(), arr.end());

    // バケットごとに桁の値を考慮してソートする
    for (int exp = 1; maxVal / exp > 0; exp *= 10) {
        std::vector<int> output(arr.size());
        std::vector<int> count(10, 0);

        for (int i = 0; i < arr.size(); i++)
            count[(arr[i] / exp) % 10]++;

        for (int i = 1; i < 10; i++)
            count[i] += count[i - 1];

        for (int i = arr.size() - 1; i >= 0; i--) {
            output[count[(arr[i] / exp) % 10] - 1] = arr[i];
            count[(arr[i] / exp) % 10]--;
        }

        arr = output;
    }
}

void stdsort(std::vector<int>& arr){
	std::sort(arr.begin(), arr.end());
}

double measureTime(std::vector<int>& arr, void (*sortingAlgorithm)(std::vector<int>&)) {
    auto start = std::chrono::high_resolution_clock::now();
    sortingAlgorithm(arr);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    return duration.count();
}

int main() {
    const int N = 10000000; // テスト用のデータサイズ
    const int numTrials = 30; // 平均を取る試行回数
    const int rndbits = __builtin_popcount(RAND_MAX);

    std::vector<int> Arr(N);

    // テスト用のランダムデータを生成
    srand(time(NULL));

    double totalFlashTime = 0.0;
    double totalRadixTime = 0.0;
    double totalstdsortTime = 0.0;

    // 平均実行時間を計測
    for (int i = 0; i < numTrials; ++i) {
        std::cout << i+1 << ": ";
        for (int i = 0; i < N; ++i) {
            int val = (rand() << rndbits) | rand();
            Arr[i] = val;
        }
        std::vector<int> flashCopy = Arr; // 元のデータをコピー
        std::vector<int> radixCopy = Arr; // 元のデータをコピー
        std::vector<int> stdsortCopy = Arr; // 元のデータをコピー

        totalFlashTime += measureTime(flashCopy, flashSort);
        totalRadixTime += measureTime(radixCopy, radixSort);
        totalstdsortTime += measureTime(stdsortCopy, stdsort);
        std::cout << "done!" << std::endl;
    }

    double avgFlashTime = totalFlashTime / numTrials;
    double avgRadixTime = totalRadixTime / numTrials;
    double avgstdsortTime = totalstdsortTime / numTrials;

    std::cout << "N = " << N << std::endl;
    std::cout << "Flash Sort の平均実行時間: \t" << avgFlashTime << "秒" << std::endl;
    std::cout << "Radix Sort の平均実行時間: \t" << avgRadixTime << "秒" << std::endl;
    std::cout << "std::sort の平均実行時間: \t" << avgstdsortTime << "秒" << std::endl;

    return 0;
}
