/*
衝突回数を記録して解析するプログラム
Stochastic method

・注意点
ts_wid（時間幅）とN_testを考慮する必要がある

σ/N_test　散乱断面積を小さくする分
N_mean*N_test　粒子数密度を増やす

・使い方
g++ N_test.cpp GaussLaguerre.cxx



*/


//ーーーーーーーーーーーーー定義ーーーーーーーーーーーーー

#include <fstream>
#include <math.h>
#include <iostream>
#include <random>
#include <cmath>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include "GaussLaguerre.h"

using namespace std;

//定数
double zeta_3 = 1.20205; // アペリーの数
double h_barc = 197.327; // 1 = ℏc ＝197.327 MeV*fm

//変化させるパラメータ
double T = 300.0; //温度 MeV
double L = 1.0/h_barc; //流体素編の一片の長さ fm=1/200(Mev)
int d_g = 16,d_q = 36; //自由度　gluon(boson):16   quark (fermion):3
int num_type_max =2;//区別する粒子の種類(0:bose、1:ferimion)

//todo N_test
double N_test = 1.0;
double ts_wid=0.1/(h_barc);//(fm)タイムステップ幅

double sigma=0.3/(pow(h_barc,2.0)*N_test);//散乱断面積(fm^2)

//ジェットの初期値
double E_jet_lab = 100.0;//(MeV)
double x_jet_ini = -0.5/h_barc;//ジェットの初期位置

double p_jet_lab_x= E_jet_lab;//massless
double v_jet_lab = p_jet_lab_x / E_jet_lab;

//実験に関して
int events=1000;//イベント回数
double seed =123;//乱数のseed
double time_max=1.0/h_barc;//(fm)時間発展の最大値

int ts_max = ceil(time_max/ts_wid);//タイムステップの数（tsループの回数使用）floorで切り捨て
int resolution = ts_max;//時間を何分割するか

//全実験の衝突回数の合計
//必要度C。全てのイベントで数える必要ないかも
int col_total_in_all_events=0;
vector<double> t_c_min_order;//tcを時間順に並べる


// todo:密度を大きくする分、ここが余分に必要
int N_spare = 10000;//粒子をイラストかするときのベクトルの要素数を余分に用意しておく

//流体素片（箱）の中の全粒子数
//{L(fm)*T(MeV) * 1 / h_barc (=1/197.327fm*MeV)} =L*T/197.327 [無次元]
double total_mean_num_boson = N_test * pow(L * T , 3.0) * d_g * zeta_3  /pow(M_PI, 2.0);//boson 
double total_mean_num_fermi = N_test *  3.0/4.0 * pow(L * T , 3.0) * d_q * zeta_3  /pow(M_PI, 2.0) ;//fermion
double N_total_mean = total_mean_num_boson + total_mean_num_fermi;
int num_par,num_bose,num_fermi;
int N_tot;

//積分した関数F(p) ※0から∞の方に変更したい。
//積分範囲 xmin ~ xmax =F(xmax) ~ F(xmin)
double xmin = 0.0001;//どこまで積分するか
double xmax = 30*T;//大きい数字に変数変換分のTをかけておく

//データを書き込むもの
std::ofstream col_times_in_1ev_out("col_times_in_1ev.dat");//t回目のイベントでの衝突回数合計

//ダブりもカウント
std::ofstream t_cols_wcounts_in_1ev_out("t_cols_wcounts_in_1ev.dat");//ダブりも含めて1evの衝突回数を出し平均自由工程と比較

//ーーーーーーーーーーーーーーーーーーーーーーーーーー

//二分法
double high_defo = xmax ;//初期値(右端)
double low_defo = 0.0;//初期値(左端)　
double cri = pow(0.1,5.0);//判定基準範囲（二分法）
double y_k;//F(P)の乱数

double d_rel_sqr=sigma/M_PI;//相対論的な粒子半径に相当。これよりbが小さかったら衝突。
double b_rel_sqr;//ローレンツ不変なのでどこからみても同じ
double t_c_CM;//重心系での最接近時間
double tc_ID;//衝突時刻の最小値のIDを求める際に使用（関数内にも使うのでグローバル変数として宣言）

//重心系
double E_tot_lab,beta_CM_x,beta_CM_y,beta_CM_z;
vector<double> beta_CM{beta_CM_x,beta_CM_y,beta_CM_z};
//Lorents transformation
vector<double> Lorentz_trans{2};//Lorentz transformation
vector<double> Lorentz_trans_2{2};//inverse transformation
vector<double> p_2_lab_bf(4),x_2_lab_bf(4);
//ーーーーーーーーーーーーー↑定義ーーーーーーーーーーーーー

//ーーーーーーーーーーーーー↓関数ーーーーーーーーーーーーー

//二分法で使用する関数
double F_integ(double p,int part_type_for_F_integ){//F(p)の関数：積分値を返す
  double x[38],xw[38];
  double sum = 0.0;
  double fx;

  GaussLaguerre::Gauss38(xmin,p,x,xw);

    if(part_type_for_F_integ == 0){

        for (int i=0;i<38;i++){

        fx = pow(x[i],2.0) / (exp(x[i]/T)-1.0);//被積分関数
        //boson

        sum += fx*xw[i];
        }
    }

    else if(part_type_for_F_integ == 1){

        for (int i=0;i<38;i++){

        fx = pow(x[i],2.0) / (exp(x[i]/T)+1.0);
        //fermion

        sum += fx*xw[i];
        }
    }
  return sum;
}
double func (double x_2,int part_type){
  
  double y;
  if(part_type == 0){//boson

    y =  F_integ(x_2,0) - y_k;

  }

  if(part_type == 1){//fermion

    y =  F_integ(x_2,1) - y_k;

  }

  return y;
}

//β(v/c)を入力して、ラピディティyを出力する関数
double beta_to_y(double beta){
    double y;
    y=(1.0/2.0)*log((1.0+beta)/(1.0-beta));
    return y;
}

//三次元の内積
double in_pro_3d(vector<double>& vec_1, vector<double>& vec_2){
 double total_3d=0;
  //空間成分
  for (int i = 0; i < 3; i++){

    total_3d += vec_1[i] * vec_2[i];
  }

  return total_3d;
}

//四次元の内積
double in_pro_4d(vector<double>& vec_1, vector<double>& vec_2){
 double space_component=0;
 double total_4d=0;
  //空間成分
  for (int i = 1; i < 4; i++) {

    space_component += vec_1[i] * vec_2[i];
  }
  //時間成分
  total_4d = (vec_1[0]*vec_2[0]) - space_component;
  return total_4d;
}

//行列の関数）　          (i行1列         ＝ i行k列 (A行列は2次元配列)　 *    k行1列の計算 (Bは1次元配列)
void cal_i_1_matrix(vector<double>& C,vector<vector<double>>& A, vector<double>& B){

    //C(n行z列) = A(n行m列) ＊　B(y行z列)を計算

    int n = A[0].size();//A行列の列数
    int m = A.size();//A行列の行数(B行列の列数と一致)

    //1,2,3, , ,i行(↓）1列目の計算
    for(int i=0; i<n; i++){

        //共通部分のkを一番内側のループにする
        //kのループで1成分分を計算(k次元の内積として、kの情報は消えて1成分になる)
        for(int k=0; k<m; k++){

            //i行は   i行k列 と  k行 の積
            C[i] += A[i][k] * B[k];

        }
    }
return;
}

//Lorents matirix reload function
vector<vector<double> > Lor_matrix{vector<double>{2},vector<double>{2}};
//Lor_matrixは関数の外で定義しておく。
void new_y_CM(double y){
Lor_matrix={
    {cosh(-y),sinh(-y)},
    {sinh(-y),cosh(-y)}
    };
    return;
}

//rotate around z-axis(transform φ)
vector<vector<double>> rot_z{vector<double>{3},vector<double>{3},vector<double>{3}};
//x-y平面の角度φを更新する関数
void new_rot_mat_phi(double phi){
rot_z={
    {cos(-phi),-sin(-phi),0},
    {sin(-phi),cos(-phi),0},
    {0,0,1}
    };
    return;
}

//rotate around y-axis(transform θ)
vector<vector<double>> rot_y{vector<double>{3},vector<double>{3},vector<double>{3}};
//z軸との角度θを更新する関数
double new_rot_mat_theta(double theta){
rot_y={
    {cos(-theta),0,sin(-theta)},
    {0,1,0},
    {-sin(-theta),0,cos(-theta)}
    };
    return 0;
}

//最小値を取り出すアルゴリズム
void min_func(vector<double>& v,double* min_tc,int* min_tc_ID){

    double min = v[0];//第0成分を暫定的な最小値に設定
    int min_ID=0;

    for (std::size_t i=0;i<v.size();i++){
        
        if(v[i]<100){//元々大きい値は弾く
            //cout << i <<"番目の候補" << v[i] << ":暫定1位:"<< min << endl;
            if (v[i] < min ) {
                min = v[i];
                min_ID = i;
            }
        }
    }
    //cout << " " << endl;
    *min_tc = min;
    *min_tc_ID = min_ID;
    return;
}

//ーーーーーーーーーーーーー↑関数ーーーーーーーーーーーーー

int main(){

double t_lab_max = time_max;
double t_step_width_scale = ts_wid;

    //cout << "ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー結果ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー" << endl;
    //cout << "" << endl;

    std::mt19937 mt(seed);//メルセンヌツイスタを使用

    //座標
    std::uniform_real_distribution<double> distr(-L/2.0,L/2.0);

    //運動量用
    std::uniform_real_distribution<double> p_dist_bose(F_integ(xmin,0),F_integ(xmax,0));//a≤f(P)≤b
    std::uniform_real_distribution<double> p_dist_fermi(F_integ(xmin,1),F_integ(xmax,1));//a≤f(P)≤b
    std::uniform_real_distribution<double> dis_theta(-1,1);//-1≤cosθ≤1
    std::uniform_real_distribution<double> dis_phi(0,2*M_PI);//0≤φ"≤"2π

cout << "mean_boson:" << total_mean_num_boson << endl;
cout << "mean_fermi:" << total_mean_num_fermi << endl;

    std::poisson_distribution<> dist_bose(total_mean_num_boson);//boson
    std::poisson_distribution<> dist_fermi(total_mean_num_fermi);//fermion

    //ーーーーー全実験を通して使う配列ーーーーー
    vector<double> col_total_vec_in_1ev(events);//実験ごとの衝突回数を入れとく箱
    vector<vector<double>> t_col_min_vec(events,vector<double>(resolution,0));//衝突時刻を入れておく箱（平均値を求める際に使用）
    //0を入れて後半のifで初期値と同じ＝値の入っていない要素は弾く
    vector<vector<double>> t_col_in_t_step(events,vector<double>(resolution));//ダブりがないかをチェック(数え漏れの数の確認)

//ダブりはいらない
int cols_total_in_all_events =0;//ダブりも含めてカウント
int cols_total_in_all_events_boson =0;//ダブりも含めてカウント
int cols_total_in_all_events_fermi =0;//ダブりも含めてカウント

//ーーーーーーーーーーーーーeventのループ開始ーーーーーーーーーーーーー

    for(int t=0; t< events;t++){

        cout << "ーーーーーーーーevents_times: " << t + 1 << " :ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー" << endl;

        //ーーーーー実験ごとに初期化する配列ーーーーーーーー
        int col_times_total =0;//t回目の実験での衝突回数
        int t_cols_wcounts_in_1ev =0;
        //(t_step→N_ID→(x,y,z))
        //[タイムステップ][ID][座標成分],粒子の座標を入れとく箱。1stepに最大1個が選抜される。
        vector<vector<vector<double>>> part_data_of_x_lab(resolution,vector<vector<double>>(N_spare,vector<double>(3)));
        
        //一回の実験での衝突位置のイラスト:
        //イラストは最悪いらない
        vector<vector<double>> col_part_data_of_x_lab(resolution,vector<double>(3));
        
        //ジェットの定義は一定じゃないので変更が必要

        // todo 初期値の入れ方
        double x_p1 = x_jet_ini;
        vector<double>p_1_lab_bf{E_jet_lab,p_jet_lab_x,0.0,0.0};//masslessでx方向にのみ移動するジェット
        vector<double>x_1_lab_bf{0.0,x_p1,0.0,0.0};

        for(int t_step_lab=0;t_step_lab<ts_max;t_step_lab++){
        
            //タイムステップごとに捨てていい配列
            //(ID,P or xyz)衝突粒子の運動量 or 座標を取り出すために必要なベクトル
            vector<vector<double>> all_p_in_1ts(0,vector<double> (4));
            vector<vector<double>> all_x_in_1ts(0,vector<double> (4));
            //ただこれは、条件分岐で衝突回数が一回でもあるなら。と使っているので絶対に消せない。
            int col_times_in_t_step=0;//tstepあたりの衝突回数をカウント
            //ダブりもカウント

            //cout << "ーーーーーーーーevents_times: " << t + 1 << " :ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー" << endl;

            //cout << "ーーーーーーーーーーーーーt_step: " << t_step_lab + 1 << " :ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー" << endl;

            //p_1が1マスすすむ：これが多分こんな単純じゃない
//ジェットは前の散乱によって運動量が変化する。
/*
前回のタイムステップの終わりにNewp1を更新しておく。
*/
x_1_lab_bf[1] = x_p1 + ( t_step_lab* t_step_width_scale * (p_1_lab_bf[1]/p_1_lab_bf[0]));

            //ーーーーーーーーーーーーー↓backgroundーーーーーーーーーーーーー

            vector<double> part_data(9);
            //{x,y,z,E,|p_x|,|p_y|,|p_z|,i,type}; i=ID 0 <= i <= N(粒子数)

            for(int type=0;type<num_type_max;type++){
                //(gluon=0,quark=1の条件分岐)ーーーーーーーーーーーーー

                if(type==0){
                    num_par = dist_bose(mt);
                    num_bose = num_par;
                    N_tot = num_bose;
                }

                else if(type==1){
                    num_par = dist_fermi(mt);
                    num_fermi = num_par;
                    N_tot = num_bose + num_fermi;
                }
                
                //gluon=0 or quark=1
                //num_parに値が入った後に0の場合
                if(num_par==0){//粒子数が0の時は衝突なし→次のタイムステップへ
                //そこで何も起こらなかったというだけなので、判定落ちと同じ扱い
                    break;//粒子割り当てのループから抜ける
                }
                
                //粒子数分tcを入れとく箱を用意(g,qをどっちも含めて衝突時刻が近い順に並べる)
                t_c_min_order.resize(N_tot);
                all_p_in_1ts.resize(N_tot,vector<double>(4));//粒子数分用意
                all_x_in_1ts.resize(N_tot,vector<double>(4));
                for(int n=0;n<num_par;n++){//粒子数分のループで割り当てとそれぞれへの判定を一気に行う

                    x_2_lab_bf[0]=0.0;//時間は0を代入
                    //x,y,zの順番に座標をランダムで振る
                    for (int k = 0; k < 3; ++k) { 

                        //割り当て用の配列は共用
                        part_data.at(k) = distr(mt);
                        x_2_lab_bf[k+1]= part_data[k];

                        //データファイル用の配列は分ける
                        if(type==0){//1周目
                            part_data_of_x_lab[t_step_lab][n][k] = part_data.at(k);                       
                            all_x_in_1ts[n][k+1] = part_data[k];
                        }

                        else if(type==1){//2周目
                            part_data_of_x_lab[t_step_lab][num_bose + n][k] = part_data.at(k);            
                            all_x_in_1ts[num_bose + n][k+1] = part_data[k];
                        }

                    }   


//消しても問題なさそう
/*
                    if(type==0){//2周目
                        all_x_in_1ev_out << h_barc*part_data_of_x_lab[t_step_lab][n][0] <<  " " << h_barc*part_data_of_x_lab[t_step_lab][n][1] << " " <<  h_barc*part_data_of_x_lab[t_step_lab][n][2] << endl;
                    }

                    if(type==1){//2周目
                        all_x_in_1ev_out << h_barc*part_data_of_x_lab[t_step_lab][num_bose +n][0] <<  " " << h_barc*part_data_of_x_lab[t_step_lab][num_bose +n][1] << " " <<  h_barc*part_data_of_x_lab[t_step_lab][num_bose +n][2] << endl;
                    }
*/

//cout << "!!!!!!!!!!          前partdataチェック用          !!!!!!!!!!:" << part_data[3]<<endl;

                    //二分法による割り当て
                    //初期化
                    double center;
                    double gap;
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

                        if (func (center,type)>=0){
                            high = center;
                        }

                        if (func (center,type)<0){
                            low = center;
                        }

                        if(func (center,type)==0||cri >= gap){

                            break;

                        }
                    }
                    //r,θ,φの割り当て                        
                    double r = center;// F(p)から二分法によってもとまったp

                    double theta = acos(dis_theta(mt));
                    double phi = dis_phi(mt) ;

                    //{x,y,z,E,|p_x|,|p_y|,|p_z|,i,type}
                    //p_x,p_y,p_zの割り当て
                    part_data.at(4) = r *sin(theta)*cos(phi);//x=rsinθcosφ
                    part_data.at(5) = r *sin(theta)*sin(phi);//y=rsinθsinφ
                    part_data.at(6) = r *cos(theta);//z=rcosθ
                    //運動量をもとにmassless粒子のエネルギーを入れる
                    part_data.at(3) = sqrt(pow(part_data.at(4),2.0)+pow(part_data.at(5),2.0)+pow(part_data.at(6),2.0));
                    part_data.at(8) = (double)type;//粒子の種類(gluon0,quark1)

                    p_2_lab_bf[0]=part_data[3];//エネルギー


 //cout << p_2_lab_bf[0] << endl;

//cout << "!!!!!!!!!!          後partdataチェック用          !!!!!!!!!!:" << part_data[3]<<endl;

                for(int i=0;i<3;i++){//初期配置と初期運動量(t=0)を決める
                    p_2_lab_bf[i+1]= part_data[i+4];//空間成分のみ                              
                }

                    if(type==0){//1周目
                        part_data.at(7) = (double)n;//粒子ID
                        all_p_in_1ts[n][1] = part_data[4];
                        all_p_in_1ts[n][2] = part_data[5];
                        all_p_in_1ts[n][3] = part_data[6];                        
                        all_p_in_1ts[n][0]=part_data[3];//エネルギー  
                    }

                    //違う型に代入するときは(型) 変数名 でキャストする

                    else if(type==1){//2周目
                        part_data.at(7) = (double)(n+ num_bose);//粒子ID
                        all_p_in_1ts[num_bose + n][1] = part_data[4];
                        all_p_in_1ts[num_bose + n][2] = part_data[5];
                        all_p_in_1ts[num_bose + n][3] = part_data[6];   
                        all_p_in_1ts[num_bose + n][0]=part_data[3];//エネルギー                    
                    }


all_x_in_1ev_out << h_barc*x_2_lab_bf[1]<< " " <<h_barc*x_2_lab_bf[2] << " " <<h_barc*x_2_lab_bf[3] << " " << 0.001*p_2_lab_bf[1]<< " " <<0.001*p_2_lab_bf[2] << " " <<0.001*p_2_lab_bf[3] <<endl;

                    //x,P,p,qを作る
                    vector<double> x_lab_bf(4),P_lab_bf(4),p_lab_bf(4),q_lab_bf(4);

                    for(int i=0;i<4;i++){

                        //相対座標
                        x_lab_bf[i] = x_1_lab_bf[i] - x_2_lab_bf[i];

                        //P=p1+p2
                        P_lab_bf[i] = p_1_lab_bf[i] + p_2_lab_bf[i];

                        //相対運動量p=p1-p2
                        p_lab_bf[i] = p_1_lab_bf[i] - p_2_lab_bf[i];

                    }

                    //内積の結果を新たな変数で置き換え
                    double x_x,x_p,x_P,p_p,P_P,p_P,q_q,x_q;
                    x_x=in_pro_4d(x_lab_bf,x_lab_bf);
                    x_p=in_pro_4d(x_lab_bf,p_lab_bf);
                    x_P=in_pro_4d(x_lab_bf,P_lab_bf);    
                    p_P=in_pro_4d(p_lab_bf,P_lab_bf);
                    P_P=in_pro_4d(P_lab_bf,P_lab_bf);
                    p_p=in_pro_4d(p_lab_bf,p_lab_bf);

                    //q,その他の計算結果の割り当て
                    for(int i=0;i<4;i++){
                    
                        q_lab_bf[i]=p_lab_bf[i]-((p_P/P_P)*P_lab_bf[i]);

                    }
                    
                    q_q=in_pro_4d(q_lab_bf,q_lab_bf);
                    x_q=in_pro_4d(x_lab_bf,q_lab_bf);

                    //新しい粒子ペアの情報を取得
                    //p1は固定。p2はbackgroundのpart_dataを使う
                    //ID番目の粒子をp_2_lab_bfとする
//それぞれの粒子のデータ
                    // cout << "=================================" << endl;
                    // cout << "粒子ID:type  " << part_data[7] << " : " << part_data[8] << endl;
                    // cout << "座標" << h_barc*part_data[0] << " " << h_barc*part_data[1] << " "<< h_barc* part_data[2] << endl;
                    // cout << "運動量"<< part_data[4] << " " << part_data[5] << " "<< part_data[6] << endl;

                b_rel_sqr=0;

                    b_rel_sqr= - x_x + (pow((x_P),2.0) / (P_P)) + (pow((x_q),2.0) / (q_q));
//確認用                    //cout << "by(P,p,q)→b_rel^2:" << pow(h_barc,2.0)*b_rel_sqr << endl;


// cout << "xx:" << x_x << endl;
// cout << "x_P:" << x_P  << endl;
// cout << "x_p:" << x_p << endl;
// cout << "x_q:" << x_q << endl;
// cout << "x:" << x_lab_bf[0] << endl;
// cout << "p:" << p_lab_bf[0] << endl;
// cout << "P:" << P_lab_bf[0] << endl;
// cout << "q:" << q_lab_bf[0] << endl;
// cout << b_rel_sqr << endl;

//todo :粒子を止めてみる


                    //散乱断面積による衝突判定
                    //if b_rel<=sqrt(散乱断面積σ-sigma/π)を満たすなら1次基準クリア
                    //最接近距離から見れば衝突しうる→あとはそれがタイムステップ内かの判定が必要

//ーーーーーーーーーーーーー衝突径数判定ーーーーーーーーーーーーー                
                    if(b_rel_sqr<=d_rel_sqr){

                        //cout << "*******b_rel_sqr<=d_rel_sqrなので衝突の可能性あり。*******" << endl;

                        //t_n_labをt_n_CMにローレンツ変換
                
                        //ローレンツ変換行列のyをy_CMに更新
                        //β(重心速度)を極座標変換(座標回転するため)
                        double beta_CM_size=0;
                        double beta_CM_size_xy=0;
                        double phi_CM=0;
                        double theta_CM=0;
                        vector<double> beta_CM_pole={beta_CM_size,phi_CM,theta_CM};

                        //CM系のパラメータ計算
                        E_tot_lab=p_1_lab_bf[0]+p_2_lab_bf[0];
                        beta_CM_x=(p_1_lab_bf[1]+p_2_lab_bf[1])/E_tot_lab;
                        beta_CM_y=(p_1_lab_bf[2]+p_2_lab_bf[2])/E_tot_lab;
                        beta_CM_z=(p_1_lab_bf[3]+p_2_lab_bf[3])/E_tot_lab;

                        double y_CM = beta_to_y(beta_CM_size);
                        new_y_CM(y_CM);

                        /*変更点：x(t_c)の第0成分からt_cを計算*/
                        //ーーーーーーーーーーーーー↓変更点ーーーーーーーーーーーーー

                        //見やすくするためにそれぞれの内積を変数でおく
                        double p1_p1,p2_p2,p1_p2,x_p1,x_p2,p2_p,p1_p,p1_P,p2_P;
                        double t_c_CM,t_1_CM,tc_t1_CM,tc_t2_CM;

                        p1_p1=in_pro_4d(p_1_lab_bf,p_1_lab_bf); 
                        p2_p2=in_pro_4d(p_2_lab_bf,p_2_lab_bf);
                        p1_p2=in_pro_4d(p_1_lab_bf,p_2_lab_bf);
                        x_p1=in_pro_4d(x_lab_bf,p_1_lab_bf);
                        x_p2=in_pro_4d(x_lab_bf,p_2_lab_bf);
                        p1_p=in_pro_4d(p_1_lab_bf,p_lab_bf);   
                        p2_p=in_pro_4d(p_2_lab_bf,p_lab_bf);
                        p1_P=in_pro_4d(p_1_lab_bf,P_lab_bf);   
                        p2_P=in_pro_4d(p_2_lab_bf,P_lab_bf);

// cout << "p1_P:" << p1_P << endl;
//  cout << "P_P:" << P_P << endl;
// cout << "x_P:" << x_P << endl;
// cout << "p2_p:" << p2_p << endl;
// cout << "x_p:" << x_p << endl;
// cout << "p2_P:" << p2_P << endl;
// cout << "p1_P:" << p1_P << endl;

//確認用
// cout << "ーーーーーーーーーーーーー" << endl;
// cout << ":" << ((p1_p*p2_P)-(p2_p*p1_P)) << endl;
// cout << ":" << (p1_P / sqrt(P_P)) * ((x_P*p2_p) - (x_p*p2_P))  << endl;

                        tc_t1_CM = (p1_P / sqrt(P_P)) * ((x_P*p2_p) - (x_p*p2_P))/((p1_p*p2_P)-(p2_p*p1_P));
//                        cout << "tc_t1_CM:" << tc_t1_CM << endl;
                        tc_t2_CM = (p2_P / sqrt(P_P)) * ((x_P*p1_p) - (x_p*p1_P))/((p1_p*p2_P)-(p2_p*p1_P));
//                         cout << "tc_t2_CM:" << tc_t2_CM << endl;

// cout << "ーーーーーーーーーーーーー" << endl;

                        // t_intvl1とtc_t2_CMは基本一致しない(軌跡との最接近時刻が同時になることはほぼない)
                        // なので別のやり方をする。

                        vector<double> x_1_lab_bf_tc(4),x_2_lab_bf_tc(4);

                        for(int i =0;i<4;i++){
//cout << "x_1_lab_bf[i]:" << i << " : " << x_1_lab_bf[i] << endl;
//cout << "p_1_lab_bf[i]:" << i << " : " << p_1_lab_bf[i] << endl;
//cout << "!!!!!!!!!!          チェック用          !!!!!!!!!!:" << tc_t1_CM<< endl;

                            x_1_lab_bf_tc[i] = x_1_lab_bf[i] + ((tc_t1_CM / (p1_P / sqrt(P_P))) * p_1_lab_bf[i]);
                            x_2_lab_bf_tc[i] = x_2_lab_bf[i] + ((tc_t2_CM / (p2_P / sqrt(P_P))) * p_2_lab_bf[i]);
//cout << "x_1_lab_bf_tc[i]:" << i << " : " << x_1_lab_bf_tc[i] << endl;

                        }

                        double tc_t1_lab = x_1_lab_bf_tc[0];
                        double tc_t2_lab = x_2_lab_bf_tc[0];

//cout << "チェック:tc_t1_lab:" << tc_t1_lab << endl;
//cout << "チェック:tc_t2_lab:" << tc_t2_lab << endl;


//                        //cout << "tc_t1_lab.tc_t2_lab:" << tc_t1_lab<< " " << tc_t2_lab<< endl;            

                        //①単純に２でわる方法
                        double tc_order_lab = (tc_t1_lab +tc_t2_lab) /2.0;
//                        //cout << "tc_order_lab:" << tc_order_lab << endl;

                        vector<double> b_rel_vec(4);

                        for(int i =0;i<4;i++){

                            b_rel_vec[i] = x_lab_bf[i] + ((tc_t1_CM / (p1_P / sqrt(P_P))) * p_1_lab_bf[i]) - ((tc_t2_CM / (p2_P / sqrt(P_P))) * p_2_lab_bf[i]);

                        }

                        b_rel_sqr = in_pro_4d(b_rel_vec,b_rel_vec);
         
//                         //cout << "***チェック***t_cから求めたb_rel^2:" << pow(h_barc,2.0)*b_rel_sqr << endl;
//                         //x_1_lab_bf_tc[0]由来のtc-t1を重心系に変換したらtc_t1_CMと一致するかの確認
//                         //cout << "tc_t1_CM:" << tc_t1_CM << endl;
//                         //cout << "tc_t2_CM:" << tc_t2_CM << endl;

                        //重心系へ移動(回転とローレンツ変換)Lor_Rot_trans.cppから抜粋
                        //回転のための極座標変換
                        beta_CM_size=sqrt(pow(beta_CM_x,2.0)+pow(beta_CM_y,2.0)+pow(beta_CM_z,2.0));
                        beta_CM_size_xy = sqrt(pow(beta_CM_x,2.0)+pow(beta_CM_y,2.0));
                        phi_CM =  acos(beta_CM_x/beta_CM_size_xy);
                        theta_CM = acos(beta_CM_z/beta_CM_size);

                        //x_1,2_lab_bf_tcの空間成分の回転
                        vector<double> x_1_lab_bf_tc_space(3),x_2_lab_bf_tc_space(3);

                        for(int j = 0;j<3;j++){

                            x_1_lab_bf_tc_space[j] = x_1_lab_bf_tc[j+1]- x_1_lab_bf[j+1];
                            x_2_lab_bf_tc_space[j] = x_2_lab_bf_tc[j+1]- x_2_lab_bf[j+1];

                        }

                        vector<double> x_1_lab_af_space(3),x_2_lab_af_space(3);

                        //回転のみbfとafで確認（結果を入れるための箱の用意）
                        //一括変換による結果を収納するベクトル
                        //衝突前)φ回転
                        vector<double> x_1_lab_bf_tc_R(3),x_2_lab_bf_tc_R(3);
                        //衝突前)φ→θ回転
                        vector<double> x_1_lab_bf_tc_RR(3),x_2_lab_bf_tc_RR(3);
                        
                        //重心系での4元運動量ベクトル(衝突前)
                        vector<double> x_1_CM_bf_tc(4),x_2_CM_bf_tc(4);//4成分で表示する

                        //θφの更新
                        new_rot_mat_phi(phi_CM);
                        new_rot_mat_theta(theta_CM);

                        //衝突前)φ回転
                        cal_i_1_matrix(x_1_lab_bf_tc_R,rot_z,x_1_lab_bf_tc_space);
                        cal_i_1_matrix(x_2_lab_bf_tc_R,rot_z,x_2_lab_bf_tc_space);

                        //衝突前)φ→θ回転 前の結果の(”C”⇦ココ,A,B)をそのまま(A,B,”C”)に代入
                        cal_i_1_matrix(x_1_lab_bf_tc_RR,rot_y,x_1_lab_bf_tc_R);
                        cal_i_1_matrix(x_2_lab_bf_tc_RR,rot_y,x_2_lab_bf_tc_R);

                        //ローレンツ変換
                        //y_CMの更新
                        new_y_CM(beta_to_y(beta_CM_size));
                        x_1_CM_bf_tc[0] = (x_1_lab_bf_tc[0]*Lor_matrix[0][0])+ (x_1_lab_bf_tc_RR[2]*Lor_matrix[0][1]);//ローレンツ変換した座標の第0成分
                        x_2_CM_bf_tc[0] = (x_2_lab_bf_tc[0]*Lor_matrix[0][0])+ (x_2_lab_bf_tc_RR[2]*Lor_matrix[0][1]);

//                        //cout << "x_1_lab_bf_tc[0]由来のtc-t1を重心系に変換したもの：：" << x_1_CM_bf_tc[0] << " " << x_2_CM_bf_tc[0] << endl;

                        //tcを並べ替える箱の用意完了

        //ーーーーーーーーーーーーー時間判定ーーーーーーーーーーーーー
                        //時刻による衝突判定
                        if((0<= tc_order_lab) && (tc_order_lab<t_step_width_scale)){//最接近時刻がタイムステップ内かどうかで判定

                            col_times_in_t_step += 1;//タイムステップ内での衝突回数のカウント
                            t_cols_wcounts_in_1ev +=1;
                            t_col_in_t_step[t][t_step_lab] += 1;//衝突圏内の粒子を全てカウント

cols_total_in_all_events+= 1;

//軌跡の設定(ダブりもカウントする)
vector<double> x_2_tc_wcounts(3);//QGP側の衝突点の座標を表す配列
for(int i=0;i<3;i++){

    x_2_tc_wcounts[i] = x_2_lab_bf[i+1] + (tc_order_lab) * (p_2_lab_bf[i+1]/p_2_lab_bf[0]);

}

                            //cout << "0 <= t_c_CM < ーーーーー****** t_step_CMなので衝突圏内。******ーーーーー" << endl;
                            //cout << " " << endl;

                            //配列のn(粒子数)番目に代入　これをN回
                            
                            if(type==0){
                                t_c_min_order[n] = tc_order_lab;//順番決める配列に代入
                                cols_total_in_all_events_boson+= 1;
                            }

                            else if(type==1){
                                t_c_min_order[num_bose + n] = tc_order_lab;//すでにgluon分が代入されているので、その次から値を入れる
                                cols_total_in_all_events_fermi+= 1;
                            }
                            //tc_ID = n;//N番目の粒子が暫定１位
                        }

                        else{//衝突時間の審査落ち

t_col_in_t_step[t][t_step_lab] +=0;

                            //順番を決める配列のn番目に最小値として選ばれないような大きい数字を入れとく
                            if(type==0){
                                t_c_min_order[n] = 100000;//順番決める配列に代入
                            }

                            if(type==1){
                                t_c_min_order[num_bose + n] = 100000;//すでにgluon分が代入されているので、その次から値を入れる
                            }
                            ////cout << "t_c_lab < 0 || t_step_lab < t_c_labなので衝突せず。" << endl;

                        }
                    }//衝突径数のif文の終わり

                    //衝突径数の審査落ち
                    else{

t_col_in_t_step[t][t_step_lab] +=0;

                        ////cout << "b_rel_sqr>d_rel_sqrなので衝突せず。" << endl;
                        if(type==0){//1周目
                            t_c_min_order[n] = 100000.0;//順番決める配列に代入
                        }

                        else if(type==1){//2周目
                            t_c_min_order[num_bose + n] = 100000.0;//すでにgluon分が代入されているので、その次から値を入れる
                        }
                        ////cout << "t_c_lab < 0 || t_step_lab < t_c_labなので衝突せず。" << endl;                   
                    }
                }//粒子数ループの終わり
            }//2種類の粒子に割り当て完了

            //cout << "-----＊＊＊＊＊num_bose:num_fermi＊＊＊＊＊------:" << num_bose << " " << num_fermi << endl;
            //cout << "合計粒子数:" << N_tot << endl;

            //ーーーーーーーーーーーーー割り当て完了ーーーーーーーーーーーーー


            

            if(N_tot == 0){//衝突可能性なし
                //cout << "ーーーーーーーーーーーーー衝突粒子なしーーーーーーーーーーーーー:" << endl;
                break;//衝突判定ループから抜けて次のタイムステップへ。
            }    
            
            //ーーーーーーーーーーーーー選ばれた粒子の変換開始ーーーーーーーーーーーーー

            double min_value;
            int min_value_ID;
            min_func(t_c_min_order,&min_value,&min_value_ID);//tcの最小値を求める関数
            //double min_value = min_func(t_c_min_order,min_value,min_value_ID);//tcの最小値を求める関数
            //これを衝突時間とする)

            //double min_value_ID = min_func_ID(t_c_min_order);//tcの最小値を求める関数
            //これを衝突時間とする)

            //cout << "min_value:" << min_value  << endl;
            //cout << "min_value_ID:" <<  min_value_ID<< endl;

//10000を入れるアルゴリズムが消えれば、後半の条件式いらない。
            if(col_times_in_t_step > 0 && min_value < 100.0){//衝突粒子が1個以上あるなら変換開始
            //1個目の条件だとtc_min=10000として出力されてしまうので、明かに大きい(100以上)の値は弾く

                //ーーー衝突点の表示ーーー
                //軌跡の設定
                //todo 多分必要ない
                vector<double> x_2_tc(3);//QGP側の衝突点の座標を表す配列
                for(int i=0;i<3;i++){
                    x_2_tc[i] = all_x_in_1ts[min_value_ID][i+1] + (min_value) * (all_p_in_1ts[min_value_ID][i+1]/all_p_in_1ts[min_value_ID][0]);
                }

                //vec_col_momentum_out << all_x_in_1ts[min_value_ID][1] << endl;
                //<<  " " << all_x_in_1ts[min_value_ID][2] <<  " " << all_x_in_1ts[min_value_ID][3] << " " <<
                //all_p_in_1ts[min_value_ID][1] <<  " " << all_p_in_1ts[min_value_ID][2] <<  " " << all_p_in_1ts[min_value_ID][3] << endl;//衝突粒子の運動量を記録

                //このタイムステップでは1回衝突が起きる。というカウント
                col_times_total += 1; 

                //cout << "t_cの最小値" << min_value << endl;
                ////cout << "粒子ID" << tc_ID << endl;
                t_col_min_vec[t][t_step_lab]=min_value;//[t]回目のイベントの時刻[t_step_lab]

                //cout << "" << endl;
                
                //ーーーーーー衝突するとみなされた2粒子の重心系へ移動ーーーーーーー
                //cout << "変換開始" << endl;










            }

            else{
                //cout << "衝突粒子なし。次のタイムステップへ" << endl;
            }

            for(int k=0;k<3;k++){
            
                //最小値の座標のデータを最接近のデータにうつす
                //イラスト化するためにここだけfmに合わせる
            ////cout << "min_ID:" << min_value_ID << endl;    
                col_part_data_of_x_lab[t_step_lab][k] = h_barc*part_data_of_x_lab[t_step_lab][min_value_ID][k];

            //cout << ":" <<  h_barc*part_data_of_x_lab[t_step_lab][min_value_ID][k] << endl;

            }

//todo 決してOK
            // x_col_data_out << col_part_data_of_x_lab[t_step_lab][0] << " " 
            // << col_part_data_of_x_lab[t_step_lab][1] << " " 
            // << col_part_data_of_x_lab[t_step_lab][2] << endl;  

            //次のタイムステップへ
        }//タイムステップループの終わり(1event終了)

        //ーーーーー1eventsが終わった後の結果の集計ーーーーーーーー
        //cout << "col_times_total:" << col_times_total << endl;
        col_total_vec_in_1ev[t] = col_times_total;

        col_times_in_1ev_out << col_times_total << endl;
       t_cols_wcounts_in_1ev_out << t_cols_wcounts_in_1ev << endl;

    }//実験回数ループの終わり

    cout << "(ダブりも含めた)全実験を合わせた衝突回数 : " << cols_total_in_all_events << endl;
    cout << "(ダブりも含めた)全実験を合わせた衝突回数 boson:fermi = " <<cols_total_in_all_events_boson << " " <<cols_total_in_all_events_fermi<< endl;
    cout << "全実験を合わせた衝突回数 : " << col_total_in_all_events << endl;

    //cout << "ーーーーー       データ確認用       ーーーーー" << endl;

    //ポアソン分布に基づいた粒子数を1つ取り出す
    //cout << "N_boson_mean:" << total_mean_num_boson << endl;
    //cout << "N_fermi_mean:" << total_mean_num_fermi << endl;
    //cout << "N_tot_mean:" << N_total_mean << endl;
    //cout << "" << endl;



    //cout << "積分値の確認:" << endl;
    //cout << "zetaから求めたF(P)(bose):" << 2 * zeta_3 * pow(T,3.0) << endl;
    //cout <<"積分値F(max):積分値F(min)" << F_integ(xmax,0) << " : "  << F_integ(xmin,0) << endl;

    //cout << "zetaから求めたF(P)(fermi):" << 3.0/2.0 * zeta_3 * pow(T,3.0) << endl;
    //cout <<"積分値F(max):積分値F(min)" << F_integ(xmax,1) << " : "  << F_integ(xmin,1) << endl;


    cout << "events:" << events << endl;
    cout << "resolution:" << resolution << endl;
    cout << "sigma:d_rel_swr  :" << sigma*pow(h_barc,2.0) << " " << d_rel_sqr*pow(h_barc,2.0)<< endl;

    double lambda = 1.0/(sigma*pow(h_barc,2.0)*(N_total_mean/pow(L*h_barc,3)));//平均自由行程(fmの次元に直す)
    double col_times_theory = 1.0 /lambda;

    cout << "平均自由行程λ(fm):" << lambda << endl;
    cout << "平均衝突回数(L/λ):" << col_times_theory << endl;
    cout << "jetの初期値(Mev):" << E_jet_lab << endl;
/*
    //グラフの表示を自動化
    FILE *gp;
    gp = popen("gnuplot","w"); //plotをopen
    fprintf(gp,"splot 'all_x_in_t_step.dat' with vectors , 'x_col_data.dat' , 'x_jet.dat' with vectors , 'vec_x_b_rel_ok.dat' with vectors\n"); 
    fprintf(gp,"pause mouse\n"); //マウス操作があるまで止まる
    fflush(gp);
    pclose(gp);
*/
}