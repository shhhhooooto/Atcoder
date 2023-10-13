//2023/09/05



#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

int main(){
double N=3;
// cin>> N;
vector<int> a(N);

//入力
// for (int i = 0; i < N; i++){

//     cin >> a[i];

// }
a[0]=2;
a[1]=7;
a[2]=4;

//ゲーム開始
int alice=0;
int bob=0;

sort(a.begin(),a.end(),greater<int>());

int loop=floor(N/2);//ループ回数
for (int i = 0; i < loop; i++){
    
    alice +=a[i];
    bob +=a[i+1];
}

if(loop%2==1){//もしNが奇数なら

    alice +=a[(int)N-1];

}

cout << alice - bob << endl;



    
}