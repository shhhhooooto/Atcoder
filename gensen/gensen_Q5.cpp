#include <iostream>
using namespace std;


int main(){

int N,A,B;
int sum=0;
cin >> N >> A >> B;

int ten_th,one_th,hun,ten,one;//それぞれの位

for(int i=1;i<N+1;i++){

    ten_th=floor(i/10000);
    one_th=floor((i-ten_th*10000)/1000);
    hun =floor((i-one_th*1000)/100);
    ten =floor((i-hun*100)/10);
    one =floor((i-ten*10)/1);

int tot = ten_th+one_th+hun+ten+one;

    if(A <= tot && tot <= B ){

        sum += i;

    }
}

cout << sum << endl;

}
