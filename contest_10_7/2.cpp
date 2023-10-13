//日付：
//目標時間：
//かかった時間：

//*ポイント

/*
*KPLTパート
K：
P：
L：
T：

vector<>をsnipetにする。
stringを比較するときは"o"じゃなくて'o'!!!
forループで入力、出力するやつをスニペット
プレイヤーと点数を対応づけるとき。どうしたらいいか
sortによって、並べ替えの情報が消えてしまった。。

*/

//ーーー　コードスタート　ーーー

#include <iostream>
#include <string>
using namespace std;

int main(){
int n;

cin >> n;

//stringがすでに複数の文字を入れられるので、1次元あればOK
vector<string> S(n);
vector<pair<int, int> > index;

//入力
for (int i = 0; i < n; ++i) {

    cin >> S[i];

}

//各行に対して勝利数をカウント
vector<int> point(n);//点数を入れておく
for (int i = 0; i < n; ++i) {

    for(int j=0;j<n;j++){

        if (S[i][j] == 'o')//勝利してるなら
        {
            /* code */
            point[i] += 1;
        }
        
    
    }

}

//勝ちが多い順番に表示

for (int i = 0; i < n; i++)
{
    /* code */
    //indexとして対応づけておく。
    index.push_back(make_pair(point[i], i));
}

sort(point.begin(), point.end(), greater<int>());


for(int j=0;j<n;j++){
    cout << index[j].second +1 << " ";
}
cout << endl;



}
