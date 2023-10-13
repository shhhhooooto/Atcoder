
#include "GaussLaguerre.h"
#include <math.h>
#include <vector>
#include <iostream>


double F_integ(double p,int part_type_for_F_integ,double T){//F(p)の関数：積分値を返す
double x[38],xw[38];
double sum = 0.0;
double fx;

double xmin = 0.0001;//どこまで積分するか
double xmax = 30*T;//大きい数字に変数変換分のTをかけておく

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

double func (double x_2,int part_type,double y_k,double T){
    
    double y;
    if(part_type == 0){//boson

        y =  F_integ(x_2,0,T) - y_k;

    }

    if(part_type == 1){//fermion

        y =  F_integ(x_2,1,T) - y_k;

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
double in_pro_3d(std::vector<double> vec_1, std::vector<double> vec_2){
    double total_3d=0;
    //空間成分
    for (int i = 0; i < 3; i++){

        total_3d += vec_1[i] * vec_2[i];
    }

    return total_3d;
}

//四次元の内積
double in_pro_4d(std::vector<double> & vec_1,std::vector<double> & vec_2){
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

//行列の関数　          i行1列         ＝ i行k列 (A行列は2次元配列)　 *    k行1列の計算 (Bは1次元配列)
void cal_i_1_matrix(std::vector<double>& C,std::vector<std::vector<double> >& A, std::vector<double>& B){

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
std::vector<std::vector<double> > Lor_matrix(2,std::vector<double>(2));
//Lor_matrixは関数の外で定義しておく。
void new_y_CM(double y){
    Lor_matrix[0][0]=cosh(-y);
    Lor_matrix[0][1]=sinh(-y);
    Lor_matrix[1][0]=sinh(-y);
    Lor_matrix[1][1]=cosh(-y);
    return;
}


//rotate around z-axis(transform φ)
std::vector<std::vector<double> > rot_z(3,std::vector<double>(3));
//x-y平面の角度φを更新する関数
void new_rot_mat_phi(double phi){
    rot_z[0][0] = cos(-phi);
    rot_z[0][1] = -sin(-phi);
    rot_z[0][2] = 0.0;
    rot_z[1][0] = sin(-phi);
    rot_z[1][1] = cos(-phi);
    rot_z[1][2] = 0.0;
    rot_z[2][0] = 0.0;
    rot_z[2][1] = 0.0;
    rot_z[2][2] = 1.0;
    return;
}

//rotate around y-axis(transform θ)
//参照渡しを使うのでベクトルで宣言
std::vector<std::vector<double> > rot_y(3,std::vector<double>(3));
//z軸との角度θを更新する関数
void new_rot_mat_theta(double theta){
    rot_y[0][0] = cos(-theta);
    rot_y[0][1] = 0.0;
    rot_y[0][2] = sin(-theta);
    rot_y[1][0] = 0.0;
    rot_y[1][1] = 1.0;
    rot_y[1][2] = 0.0;
    rot_y[2][0] = -sin(-theta);
    rot_y[2][1] = 0.0;
    rot_y[2][2] = cos(-theta);
    return;
}






//ーーーーー極⇄直行座標変換の関数ーーーーー
//直交座標のベクトルC→極座標Pベクトルにする関数
double Cart_to_Pole(std::vector<double>& P,std::vector<double>& C){
std::vector<double> C_space(3);//Cの空間成分ベクトルを用意 
//4成分の行列のとき、空間成分のみを変換
if(C.size()==4){
    for(int i=0;i<3;i++){

//Cの空間成分(1,2,3)をnewCの各成分に代入(0,1,2) 
        C_space[i]=C[i+1];
    }
}

else if(C.size()==3){
    for(int i=0;i<3;i++){
        C_space[i]=C[i];//元々空間成分のベクトルならそのまま代入
    }
}

double C_size,C_size_xy,phi_C,theta_C;
//P(|C|,θ_C,φ_C)
C_size =sqrt(pow(C_space[0],2)+pow(C_space[1],2)+pow(C_space[2],2));
C_size_xy = sqrt(pow(C_space[0],2)+pow(C_space[1],2));
theta_C = acos(C_space[2]/C_size);
double sign;//符号関数
    if(C_space[1]>=0.0){
        sign=1.0;
    }
    else if(C_space[1]<0.0){
        sign=-1.0;
    }

phi_C = sign*acos(C_space[0]/C_size_xy);

P[0]=C_size;
P[1]=theta_C;
P[2]=phi_C;

    return 0;
}

//極座標のベクトルP→直交座標のベクトルCに変換する
double Pole_to_Cart(std::vector<double>&C ,std::vector<double>& P){

std::vector<double> C_space(3);//代入先のベクトルの空間成分ベクトルを用意 
//P(|C|,θ,φ)
C_space[0]=P[0]*sin(P[1])*cos(P[2]);
C_space[1]=P[0]*sin(P[1])*sin(P[2]);
C_space[2]=P[0]*cos(P[1]);

//Cが4成分のとき、変換した空間成分を代入
if(C.size()==4){//4元ベクトルのとき
    for(int i=0;i<3;i++){
    //変換した空間成分をCの空間成分へ代入
        C[i+1]=C_space[i];
    }
}

else if(C.size()==3){//3元ベクトルのとき
    for(int i=0;i<3;i++){
        C[i]=C_space[i];
    }
}

    return 0;

}


//重心系へと変換させる関数
//変換前を受け取って中で変換後に代入する
void trans_to_CM(std::vector<double> & p_j_l_b, std::vector<double> p_Q_l_b,
    std::vector<double> & p1_b_CM,std::vector<double> & p2_b_CM,std::vector<double> & beta_CM_pole){

//ーーー重心系の情報を取得ーーー
    //重心系関係のパラメータを初期化
    double E_tot_lab,beta_CM_size,phi_CM,theta_CM;
    std::vector<double> beta_CM_cart(3);//betaの計算で使用

    //CM系(直交座標ver)のパラメータ計算
    E_tot_lab=p_j_l_b[0]+p_Q_l_b[0];
    for (int j = 0; j < 3; j++){
        //空間成分(j+1)を移し替える。
        beta_CM_cart[j] = (p_j_l_b[j+1]+p_Q_l_b[j+1])/E_tot_lab;
    }

    Cart_to_Pole(beta_CM_pole,beta_CM_cart);
    beta_CM_size=beta_CM_pole[0];
    theta_CM=beta_CM_pole[1];
    phi_CM=beta_CM_pole[2];

    //0-0時間成分と空間成分のベクトルを分ける
    //エネルギーは回転に関与しないので、空間成分(3成分)のみ用意
    std::vector<double> p_jet_lab_bf_space(3),p_QGP_lab_bf_space(3);

    //衝突前)1-1 φ回転
    std::vector<double> p_1_bf_R_space(3),p_2_bf_R_space(3);
    //衝突前)1-2 φ→θ回転
    std::vector<double> p_1_bf_RR_space(3),p_2_bf_RR_space(3);

    //0-0 4成分を3成分へ
    for(int j = 0;j<3;j++){
    //空間の0成分に4元ベクトルの1成分を入れる
        p_jet_lab_bf_space[j] = p_j_l_b[j+1];
        p_QGP_lab_bf_space[j] = p_Q_l_b[j+1];
    }

    //衝突前)1-1 φの更新と回転
    new_rot_mat_phi(phi_CM);
    cal_i_1_matrix(p_1_bf_R_space,rot_z,p_jet_lab_bf_space);
    cal_i_1_matrix(p_2_bf_R_space,rot_z,p_QGP_lab_bf_space);

    //衝突前)1-2 θの更新とφ→θ回転 前の結果の(”C”⇦ココ,A,B)をそのまま(A,B,”C”)に代入
    new_rot_mat_theta(theta_CM);
    cal_i_1_matrix(p_1_bf_RR_space,rot_y,p_1_bf_R_space);
    cal_i_1_matrix(p_2_bf_RR_space,rot_y,p_2_bf_R_space);

    //cout << "φ→θ回転後p1:" << p_1_bf_RR_space[0] << " " << p_1_bf_RR_space[1] << " " << p_1_bf_RR_space[2] << endl;
    //cout << "φ→θ回転後p2:" << p_2_bf_RR_space[0] << " " << p_2_bf_RR_space[1] << " " << p_2_bf_RR_space[2] << endl;

    //ーーーーーー2-1.2-2 重心系への変換ーーーーーー

    //y_CMの更新 (重心系のラピディティに対してローレンツ変換するため)
    new_y_CM(beta_to_y(beta_CM_size));
    //第1,2成分は変換させずにそのまま代入
    //左のベクトルは4成分.成分に注意
    p1_b_CM[1]=p_1_bf_RR_space[0];
    p1_b_CM[2]=p_1_bf_RR_space[1];
    p2_b_CM[1]=p_2_bf_RR_space[0];
    p2_b_CM[2]=p_2_bf_RR_space[1];
    //p1
    p1_b_CM[0]= (p_j_l_b[0]*Lor_matrix[0][0])+ (p_1_bf_RR_space[2]*Lor_matrix[0][1]);
    p1_b_CM[3]= (p_j_l_b[0]*Lor_matrix[1][0])+ (p_1_bf_RR_space[2]*Lor_matrix[1][1]);
    //p2
    p2_b_CM[0]= (p_Q_l_b[0]*Lor_matrix[0][0])+ (p_2_bf_RR_space[2]*Lor_matrix[0][1]);
    p2_b_CM[3]= (p_Q_l_b[0]*Lor_matrix[1][0])+ (p_2_bf_RR_space[2]*Lor_matrix[1][1]);
    //cout << "Lorentz変換後p1(E1,p1z):" << p_1_bf_CM[0] << "," << p_1_bf_CM[3] << endl;
    //cout << "Lorentz変換後p2(E2,p2z):" << p_2_bf_CM[0] << "," << p_2_bf_CM[3] << endl;
    return;
}

void reverse_trans_to_CM(std::vector<double> & p_1_af_CM, std::vector<double> p_2_af_CM,
    std::vector<double> & p_jet_lab_af,std::vector<double> & p_QGP_lab_af,std::vector<double> & beta_CM_pole){

    //衝突後)3-1 逆ローレンツ変換後
    //3-2 このタイミングで4成分を3成分へ
    std::vector<double> p_1_af_RR_space(3),p_2_af_RR_space(3);

    //衝突後)4-1 φ→θ→　col　→-θ回転
    std::vector<double> p_1_af_R_space(3),p_2_af_R_space(3);
    //衝突後)4-2 φ→θ→ col　→-θ→-φ回転 //空間成分のみ回転行列をかける
    std::vector<double> p_1_af_space(3),p_2_af_space(3);

    //y_CMの更新(-betaで逆ローレンツ変換。実験系に戻る)
    new_y_CM(beta_to_y(-beta_CM_pole[0]));

    //第1,2成分は変換させずにそのまま代入。
    p_1_af_RR_space[0]=p_1_af_CM[1];
    p_1_af_RR_space[1]=p_1_af_CM[2];
    p_2_af_RR_space[0]=p_2_af_CM[1];
    p_2_af_RR_space[1]=p_2_af_CM[2];


    //逆ローレンツ変換
    //エネルギーは回転に関与しないので、直接衝突後の実験系へ
    //回転する空間成分のみ回転用のベクトルにうつす。
    p_jet_lab_af[0]= (p_1_af_CM[0]*Lor_matrix[0][0])+ (p_1_af_CM[3]*Lor_matrix[0][1]);
    p_1_af_RR_space[2]= (p_1_af_CM[0]*Lor_matrix[1][0])+ (p_1_af_CM[3]*Lor_matrix[1][1]);

    p_QGP_lab_af[0]= (p_2_af_CM[0]*Lor_matrix[0][0])+ (p_2_af_CM[3]*Lor_matrix[0][1]);
    p_2_af_RR_space[2]= (p_2_af_CM[0]*Lor_matrix[1][0])+ (p_2_af_CM[3]*Lor_matrix[1][1]);

    //θφの更新(逆回転なので-をつけてφ→θの順番に入れ替える)
    new_rot_mat_theta(-beta_CM_pole[1]);
    new_rot_mat_phi(-beta_CM_pole[2]);


    //衝突後)φ→θ→　col　→ç回転
    cal_i_1_matrix(p_1_af_R_space,rot_y,p_1_af_RR_space);
    cal_i_1_matrix(p_2_af_R_space,rot_y,p_2_af_RR_space);

    //衝突後)φ→θ→ col　→-θ→-φ回転
    cal_i_1_matrix(p_1_af_space,rot_z,p_1_af_R_space);
    cal_i_1_matrix(p_2_af_space,rot_z,p_2_af_R_space);

    for(int j = 0;j<3;j++){
    //4元ベクトルの第1成分に空間の0成分を入れる
        p_jet_lab_af[j+1] = p_1_af_space[j];
        p_QGP_lab_af[j+1] = p_2_af_space[j];
    }

    return;
}

void conserve_check_func(std::vector<double> & p_jet_lab_bf,std::vector<double> & p_QGP_lab_bf,
    std::vector<double> & p_jet_lab_af,std::vector<double> & p_QGP_lab_af){

    //std::cout << "-- 保存量の確認 --" << std::endl;
    std::vector<double> conserved_bf(4);//変換前の保存量P=(E,Px,Py,Pz)
    std::vector<double> conserved_af(4);//変換後

    //std::cout << "前)P=(E,Px,Py,Pz)" << std::endl;
    // for (int j = 0; j < 4; ++j) {
    //     conserved_bf[j] =  p_jet_lab_bf[j] + p_QGP_lab_bf[j];
    //     std::cout << conserved_bf[j] << "; " ;
    // }
    //std::cout << std::endl;

    //std::cout << "" << std::endl;
    //std::cout << "↓↓↓衝突後" << std::endl;

    //std::cout << "後)P=(E,Px,Py,Pz)" << std::endl;
    //for (int j = 0; j < 4; ++j) {
        //conserved_af[j] =  p_jet_lab_af[j] + p_QGP_lab_af[j];
        //std::cout << conserved_af[j] << "; " ;
    //}
    //std::cout << std::endl;

    //保存量の確認
    // int conserved_comp=0;//その成分が保存していたらカウント
    // for (int j = 0; j < 4; j++)
    // {
    //     //完全に一致はしないので誤差基準を設定
    //     if((conserved_bf[j] - conserved_af[j]) < pow(0.1,5) ){
    //         conserved_comp += 1;
    //     }
    // }
    // //4成分とも保存していたら
    // if(conserved_comp == 4){
    //     //check_all_col += 1;//保存○の衝突にカウント
    //     //check_1ev += 1;//保存○のイベントにカウント
    // }



    return ;
}