
//設計図

//todo 変更点
/*
セルあたりの平均値にする

todo 注意点
・多分、粒子をセル内に巻くのはPによる判定後でいい。逆にPの前に巻いても意味がない
・0≦P≦1となるようにパラメータを決める

todo パラメータの決め方について→小さくするときは慎重に。型の違いなども影響してしまう？
x_wid→0をしすぎると、cell内の平均粒子が極端に小さくなり、小数点の計算でコンピュータの影響が出そう
x_widのスケールをa変えたら、N_testをa^3で変える。そうすれば、cell内の粒子密度を一定に保てる

ts_wid→0をしすぎると、P_collが極端に小さくなり、0~1の乱数のアタリが悪くなりそう。

todo 考えること
・セルをjetから前後にとるか、前だけにとるか

todo 間違えやすいところ
単位換算→h_barc

todo 確認
x_wid=1fmにして、N/1fm^3=n_boxの平均粒子数と一致した→サンプル自体は機能している
→これを細かくするとどうなる？
→かなりうまくいっている。理論値と一致している
*サンプリングには異常なし！

todo 懸念てん
v_relの取り方

todo 論文
フレームによる影響など、かなり詳しく載っている
https://arxiv.org/pdf/hep-ph/0406278.pdf


使い方
g++ Stochastic.cpp GaussLaguerre.cxx Function_for_col_no_cout.cpp

*/

//*include
#include <fstream>
#include <math.h>
#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include "GaussLaguerre.h"
#include "Function_for_col.h"

using namespace std;

//定数
double zeta_3 = 1.20205; // アペリーの数≒ζ(3)
double h_barc = 197.327; // 1 = ℏc ＝197.327MeV*fm

//*パラメータ
int events =1000;
double seed = 123;

//todo 注意するパラメーター
//! Pは必ず1より小さい必要がある
double N_test = 1000000.0; //こいつでσを小さくしてあげないと、セルよりも大きくなってしまう。
double ts_wid = 0.0001/h_barc;//(fm/c)
double x_wid = 0.001/h_barc;//(fm) →0.01くらいが適正？10/14 0.001にするとN_testを1000000くらいにする必要がある
double v_cell = pow(x_wid,3); 
double time_max=1.0/h_barc;//(fm)時間発展の最大値
int ts_max = floor(time_max/ts_wid);//タイムステップの数（tsループの回数使用）floorで切り捨て

double T = 300.0; //温度 MeV
double L = 1.0/h_barc; //流体素編の一片の長さ fm=1/200(Mev)
int d_g = 16,d_q = 36; //自由度　gluon(boson):16   quark (fermion):36
int num_type_max =2;//区別する粒子の種類(0:bose、1:ferimion)
double sigma=0.3/(N_test*pow(h_barc,2.0));//散乱断面積(fm^2)

//ジェットの初期値
double E_jet_lab = 100.0;//(MeV)
double p_jet_lab_x= E_jet_lab;
double v_jet_lab = p_jet_lab_x / E_jet_lab;

//流体素片（箱）の中の全粒子数
//{L(fm)*T(MeV) * 1 / h_barc (=1/197.327fm*MeV)} =L*T/197.327 [無次元]
//todo cellあたりに変更する
double mean_n_bose_cell = N_test* pow(x_wid * T , 3.0) * d_g * zeta_3  /pow(M_PI, 2.0) ;//boson
double mean_n_fermi_cell = N_test* 3.0/4.0 * pow(x_wid * T , 3.0) * d_q * zeta_3  /pow(M_PI, 2.0) ;//fermion
double mean_n_tot_cell = mean_n_bose_cell +mean_n_fermi_cell;

//全イベントを通して使う配列・変数
double P_coll_tot_for_mean = 0;//最後に判定回数で割って、P_collの平均値を出す
double v_sweep_tot_for_mean = 0;//平均値を出すため
int num_coll_judge_all_ev =0;//判定回数
int sampled_n_tot_all =0;//サンプルされた粒子数
int coll_times_all_ev_single =0;
int coll_times_all_ev_wcoutns =0;

//積分範囲 xmin ~ xmax =F(xmax) ~ F(xmin)
double xmin = 0.0001;//どこまで積分するか
double xmax = 30*T;//大きい数字に変数変換分のTをかけておく

//二分法
double high_defo = xmax ;//初期値(右端)
double low_defo = 0.0;//初期値(左端)　
double cri = pow(0.1,5.0);//判定基準範囲（二分法）
double y_k;//F(P)の乱数

//一様分布
std::mt19937 mt(seed);//メルセンヌツイスタを使用
std::uniform_real_distribution<double> p_dist_bose;
std::uniform_real_distribution<double> p_dist_fermi;
std::uniform_real_distribution<double> dis_theta(-1,1);//-1≤cosθ≤1
std::uniform_real_distribution<double> dis_phi(0,2*M_PI);//0≤φ"≤"2π
// todo x_widに変更
std::uniform_real_distribution<double> x_dist(-x_wid/2.0,x_wid/2.0);//座標 
std::uniform_real_distribution<double> P_coll_dist(0,1.0);//P_collに使用

//ポアソン分布
std::poisson_distribution<> dist_bose_cell(mean_n_bose_cell);//boson粒子数
std::poisson_distribution<> dist_fermi_cell(mean_n_fermi_cell);//fermion粒子数

//todo p,xの割り当て関数
void p_assign(int type,vector<double>&p_QGP_cell){

    //二分法による割り当て
    //初期化
    double center,gap;
    double high = high_defo;
    double low = low_defo;

    //運動量の大きさの決定
    if(type==0){
        y_k = p_dist_bose(mt);
    }

    else if(type==1){
        y_k = p_dist_fermi(mt);   
    }

    for(int i=0;i<30;i++){ //二分法で乱数(y_k)に対応する解Pを求める

        gap= abs(high-low);
        center = (high + low)/2;

        if (func (center,type,y_k,T)>=0){
            high = center;
        }

        if (func (center,type,y_k,T)<0){
            low = center;
        }

        if(func (center,type,y_k,T)==0||cri >= gap){
            break;
        }
    }
    //r,θ,φの割り当て 
    vector<double> p_QGP_pole(3);                       
    p_QGP_pole[0] = center;// F(p)から二分法によってもとまったp
    p_QGP_pole[1] = acos(dis_theta(mt));
    p_QGP_pole[2] = dis_phi(mt);
    cout << "|p_QGP|:" << p_QGP_pole[0] << endl;
    Pole_to_Cart(p_QGP_cell,p_QGP_pole);
    p_QGP_cell[0] = center;//masslessなので、|p|=E
}

void x_assign(vector<double>&x_QGP_cell){

    for(int i=0;i<3;i++){
    
    x_QGP_cell[i] = x_dist(mt);

    }

}

//散乱が起こる(1)か起こらない(0)かを判定
bool coll_judge(vector<vector<double> >&p_jet_1ev,int ts_lab,vector<vector<double> >&p_cell, int index){
    vector<double> P(4);//運動量の和
    for(int i=0;i<4;i++){
        P[i] = p_jet_1ev[ts_lab][i]+p_cell[index][i];
    }
    double s = in_pro_4d(P,P);
    double v_rel = s/(2*(p_jet_1ev[ts_lab][0]*p_cell[index][0]));
    double v_sweep = v_rel*ts_wid*sigma;//断面積がsweepする体積
    v_sweep_tot_for_mean += v_sweep;

    cout << "s:" << s << endl;
    cout << "v_rel:" << v_rel << endl;
    cout << "v_sweep:" << v_sweep << endl;

    cout << "dx^3:" << v_cell << endl;

    //セルよりもv_rel*sigmaは小さい必要がある→N_test→∞,ts_wid→0,x→大
    if(v_cell< v_rel*ts_wid*sigma){
        cout << "!!!!!!!!!!          描像が破綻している          !!!!!!!!!!:" << endl;
    }

    double P_coll = (v_sweep)/(v_cell);
    P_coll_tot_for_mean += P_coll;
    bool coll_judge=0;
    double random = P_coll_dist(mt);
    cout << "P_coll:" << P_coll << endl;
     cout << "random:" << random << endl;
    
    //この範囲に入れば衝突が起きる
    if(P_coll >= random){
        coll_judge = 1;
    }
    return coll_judge;
}   

int main(){




//*ーーーーーーーーーーーーー乱数系ーーーーーーーーーーーーー
//分布関数の範囲に積分関数を使用しているのでmain関数の中で定義する必要がある。
p_dist_bose.param(std::uniform_real_distribution<>::param_type(F_integ(xmin,0,T),F_integ(xmax,0,T)));
//std::uniform_real_distribution<double> p_dist_fermi(F_integ(xmin,1),F_integ(xmax,1));//運動量P a≤f(P)≤b
p_dist_fermi.param(std::uniform_real_distribution<>::param_type(F_integ(xmin,1,T),F_integ(xmax,1,T)));
//?　これって積分する必要あるんだっけ？普通に解析会でよくない？それならどこでも宣言できる


//todo t→eに変更
for(int e=0;e<events;e++){
cout << "ーーーーーーーーーー event =" << e<< " ーーーーーーーーーー" << endl;

//イベントごとに初期化する配列・変数
vector<vector<double> > x_jet_1ev(ts_max+2,vector<double>(4));
vector<vector<double> > p_jet_1ev(ts_max+2,vector<double>(4));
int coll_times_1ev_single_count=0;
int coll_times_1ev_wcounts=0;
int sampled_n_bose_tot_1ev =0;
int sampled_n_fermi_tot_1ev =0;


//ーーーーーーーーーーーーー流体素片に注入するジェットの初期値ーーーーーーーーーーーーー
//時刻0(初期値)における1つ目(ID=0)のジェットの情報
//ージェットの初期座標ー
x_jet_1ev[0][0]= 0.0;//時間成分 fm
x_jet_1ev[0][1]= 0.0;//x fm
x_jet_1ev[0][2]= 0.0;//y
x_jet_1ev[0][3]= 0.0;//z
//ージェットの初期運動量ー
//*(massless)
p_jet_1ev[0][0]= E_jet_lab;//エネルギー MeV
p_jet_1ev[0][1]= p_jet_lab_x;//px MeV
p_jet_1ev[0][2]= 0.0;//py
p_jet_1ev[0][3]= 0.0;//pz

    for(int ts_lab=0;ts_lab<ts_max;ts_lab++){
    cout << "ーーーーー ts =" << ts_lab << " ーーーーー" << endl;
        //tsごとに初期化する配列・変数
        bool coll_happen=0;//collが起こるなら1を代入。起こるか起こらないかが知りたい。

        //todo jetの初期座標を元にセルを考える
        //もはや考えなくてもいい？セルの大きさは固定だから毎回、1jetにつき固定の1セルを考えればいい？
        //V_cellは固定だなこりゃ
        

        //todo セルの中に粒子をサンプルする(cellあたり)
        //tsごとに初期化
        double n_bose_cell = dist_bose_cell(mt);
        double n_fermi_cell = dist_fermi_cell(mt);

        //もしセル内に粒子があるなら
        int n_cell = n_bose_cell + n_fermi_cell;
        sampled_n_bose_tot_1ev += n_bose_cell;
        sampled_n_fermi_tot_1ev += n_fermi_cell;

        // cout << "n_bose_cell:" << n_bose_cell << endl;
        // cout << "n_fermi_cell:" << n_fermi_cell << endl;
        
        if(n_cell>0){//ただし、本当はここが1になるのが理想
            
            vector<vector<double> > p_bose_cell(n_bose_cell,vector<double>(4));//bose用の配列を作る
            vector<vector<double> > p_fermi_cell(n_fermi_cell,vector<double>(4));//bose用の配列を作る
            //cell中の粒子数分ループ
            for(int i=0;i<n_bose_cell;i++){
                num_coll_judge_all_ev+=1;

                //todo pの割り当て
                vector<double> p_bose_assign(4);//割り当てにのみ使用
                p_assign(0,p_bose_assign);//i番目のボーズ粒子に割り当てる
                for(int k=0;k<4;k++){
                    p_bose_cell[i][k]= p_bose_assign[k];//空間成分に代入
                }
                
                if(coll_judge(p_jet_1ev,ts_lab,p_bose_cell,i)){//衝突が起こるなら
                    coll_times_1ev_wcounts +=1;
                    coll_happen =1;
                }
            }

            for(int i=0;i<n_fermi_cell;i++){
                num_coll_judge_all_ev+=1;

                vector<double> p_fermi_assign(4);//割り当てにのみ使用
                p_assign(1,p_fermi_assign);//i番目のボーズ粒子に割り当てる
                for(int k=0;k<4;k++){
                    p_fermi_cell[i][k]= p_fermi_assign[k];
                }

                if(coll_judge(p_jet_1ev,ts_lab,p_fermi_cell,i)){//衝突が起こるなら
                    coll_times_1ev_wcounts +=1;
                    coll_happen =1;
                }
            }
        }

        //セル内に粒子がいないなら
        else{
        
        

        }

        if(coll_happen){//衝突が起きるtsなら
            //single_countにカウント(ダブりをカウントしない)
            coll_times_1ev_single_count +=1;
        }

        //todo 軌跡を進める（今回は散乱なし）
        //この場合はx方向のみの移動でOK
        //次のタイムステップに代入
        //座標と運動量を移し替える
        for(int i=0;i<3;i++){
            x_jet_1ev[ts_lab+1][i+1] = x_jet_1ev[ts_lab][i+1] + (ts_wid * (p_jet_1ev[ts_lab][i+1]/p_jet_1ev[ts_lab][0]));        
            p_jet_1ev[ts_lab+1][i+1] = p_jet_1ev[ts_lab][i+1];
        }
        p_jet_1ev[ts_lab+1][0] = p_jet_1ev[ts_lab][0];

        cout << "x_jet:" << x_jet_1ev[ts_lab][1]*h_barc << endl;


    }
    //全イベントの集計に保存
    coll_times_all_ev_wcoutns += coll_times_1ev_wcounts;
    coll_times_all_ev_single += coll_times_1ev_single_count;
    sampled_n_tot_all += sampled_n_bose_tot_1ev + sampled_n_fermi_tot_1ev;

    //todoーーーーーー1イベント終了時のデータ集計
    cout << "ts_max:" << ts_max << endl;
    cout << "" << endl;
    cout << "---       衝突回数(1ev)       ---" << endl;
    cout << "シミュレーション結果" << endl;
    cout << "coll_times_single_1ev:" << coll_times_1ev_single_count << endl;
    cout << "coll_times_w_1ev:" << coll_times_1ev_wcounts << endl;
    cout << "理論値(L/λ):" << L*sigma*pow(h_barc*(L/x_wid),3)*mean_n_tot_cell << endl;

    cout << "" << endl;
    cout << "sampled_bose_1ev_tot:" << sampled_n_bose_tot_1ev << " 理論値:" << mean_n_bose_cell*ts_max << endl;
    cout << "sampled_fermi_1ev_tot:" << sampled_n_fermi_tot_1ev << " 理論値:" << mean_n_fermi_cell*ts_max << endl;
    cout << "" << endl;

}

//全イベント終了時のデータ集計
cout << "mean_n_boson_cell:" << mean_n_bose_cell << endl;
cout << "mean_n_fermi_cell:" << mean_n_fermi_cell << endl;
cout << "sigma:" << sigma << endl;

cout << "" << endl;
cout << "mean_N_boson_box(N_testを考えないとき):" << mean_n_bose_cell*pow((L/x_wid), 3)/N_test << endl;
cout << "mean_N_fermi_box(N_testを考えないとき):" << mean_n_fermi_cell*pow((L/x_wid), 3)/N_test << endl;
cout << "" << endl;

cout << "v_sweep_mean:" << v_sweep_tot_for_mean/num_coll_judge_all_ev << endl;
cout << "v_cell:" << v_cell << endl;
cout << "P_collの理想値(v_sweep_mean/v_cell):" << v_sweep_tot_for_mean/(num_coll_judge_all_ev*v_cell) << " →これが1以上ならアウト:" << endl;
if(v_sweep_tot_for_mean/(num_coll_judge_all_ev*v_cell)>1){
    cout << "!!!!!!!!!!                    !!!!!!!!!!:" << endl;
    cout << "!!!!!!!!!!        Pが1を超えている         !!!!!!!!!!:" << endl;
    cout << "!!!!!!!!!!                    !!!!!!!!!!:" << endl;
}
cout << "P_coll_mean:" << P_coll_tot_for_mean/num_coll_judge_all_ev << endl;


cout << "cell内の平均粒子数(結果)=サンプルした粒子の合計/サンプルした数(タイムステップの数*イベント数):" << (double)sampled_n_tot_all/((double)ts_max*(double)events) << " = cell内の平均粒子数(理論値):" <<mean_n_tot_cell << endl;


cout << "" << endl;
//w,sの衝突回数とサンプル数
cout << "1.サンプルした粒子数の合計:" << sampled_n_tot_all << endl;
cout << "2.セル内の衝突確率平均(P_coll_mean):" << P_coll_tot_for_mean/num_coll_judge_all_ev << endl;
cout << "衝突回数の理想値(1.*2.=サンプルされたうち,衝突と判定されるもの=衝突回数):   →   " << sampled_n_tot_all* P_coll_tot_for_mean/num_coll_judge_all_ev<< endl;

cout << "---     衝突回数(all_ev)     ---" << endl;
cout << "coll_times_all_ev_single:   →   " << coll_times_all_ev_single << endl;
cout << "coll_times_all_ev_wcounts:   →   " << coll_times_all_ev_wcoutns << endl;
cout << "ダブりが出た回数:" << coll_times_all_ev_wcoutns - coll_times_all_ev_single << endl;

}
