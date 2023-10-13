//2023/09/07

#include <iostream>
#include <algorithm>
using namespace std;

int main(){
int N;
cin >> N;


int d[N];
int max_pile=1;//重ねられる最大値



for(int i=0;i<N;i++){

    cin >> d[i];

}

//sort()して
sort(d, d+N,greater<int>());

//a>bならOK

for (int i = 0; i < N; i++){
    if(i != N-1){
        if (d[i+1]<d[i]){
            
            max_pile +=1;        
            // cout << "i:max:d:" << i << ":" << max_pile <<":" << d[i] << ":" << d[i+1] << endl;
        }

        else if (d[i+1]==d[i]){
            //同じなら何もしない
            // cout << "Hello" << endl;
        }
    }
    
    else{
        // cout << "end:" << i <<":"<< d[i] << endl;
        break;

    }
    

}



cout << max_pile << endl;
    
}
