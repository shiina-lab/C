/*****************************************************************
 * 
 * MRR L1 CDFファイルから降水粒子群の落下速度と時間毎のドップラースペクトルを取得し
 * 最小二乗法を用いて降水粒子の発生高度を予測する
 *
 * ===============================================================
 * メモ
 * ---------------------------------------------------------------
******************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"CDF_L1_MRR.h"
#include"CDF_L2_MRR.h"

#define CDF_DATA_LEVEL         "1"
#define CDF_DATA_VERSION      "02"
#define CDF_DATA_SUFFIX        ""

#define PATH_MAX               (1024)

//#define HEIGHT                  35 //[m]
#define DH1                      "100" //[m]
#define DH2                      "050" //[m]
#define DH3                      "200" //[m]
#define DT                      10 //[sec]
#define FNUM                    16

//h0, t0
#define H0      2082.4244         //発生高度
#define T0      28.295033         //発生時間

#define V_max  12.0
#define V_min   2.0

//** 1ヶ月間のみ対応*********//

#define DIR     "/Users/Shiina/Desktop/cdf-workspace2017/CDF-src/mrr_tani/"
#define DATE     "0825"
//#define DATE     "test"
//variable
//how long
#define DAYS                    1 //[日]
#define HOURS                   1  //[hour]
#define MINUTES                 60  //[minute]

#define TAU                     144 //range point

//prototype declaration
double  fx(double tau, double ha, double h0, double t0);

double  Zc(double tau, double vel, double ha, double h0, double dv);




int main(int argc,const char *argv[]){
    //when
    //start
    float YEAR_START         =     2017.0,
          MONTH_START        =     8.0,
          DATE_START         =     25.0,
          HOUR_START         =     10.0,
          MINUTE_START       =     0.0,
          SECOND_START=0,
    //end
          YEAR_END          =      2017.0,
          MONTH_END         =      8.0,
          DATE_END          =      25.0,
          HOUR_END          =      12.0,           // max +6 hours stack size
          MINUTE_END        =      0.0,
          SECOND_END=0;
    //*DAYS or HOURS******************
    //int data_sum = 60 * 60 * 24 * DAYS /10;     //DAYS
    //int data_sum = 60 * 60 * HOURS /10;           //HOURS
    int data_sum = ( 60*(60*HOUR_END+MINUTE_END)-60*(60*HOUR_START+MINUTE_START) ) /10;
    //********************************
    
    double   se=0;
////////////////////////////////////////
    
  // ************************************
  // MRR データ
  // ------------------------------------
  CDF_L1_MRR__DATA mrr1,mrr2,mrr3;
  // ************************************
  // 作業用変数
  // ------------------------------------
  CDFid            id1 = NULL;
  CDFid            id2 = NULL;
  CDFid            id3 = NULL;
  
    char             dirname1 [PATH_MAX];
    char             dirname2 [PATH_MAX];
    char             dirname3 [PATH_MAX];

  char             filename1[PATH_MAX];
  char             filename2[PATH_MAX];
  char             filename3[PATH_MAX];
  int              recNum;
  const char      *path;
  CDFstatus        status;
  double           year,month,day,hour,minute,second;
  int              YY,MM,DD,hh,mm,ss;
  int              Y_,M_,D_;
  
    double  sec, min, hou;
    double h,m;
    int count=0;//ループ回数をカウント
    int i,j,k,l;
    FILE *fp,*gp;
    //**************************************************************************************************************
    CDF__TIME_TT2000 epoch_start    = CDF_TT2000_from_UTC_parts(YEAR_START, MONTH_START, DATE_START, HOUR_START, MINUTE_START, 0.0,  0.0,  0.0,  0.0);
    CDF__TIME_TT2000 epoch_end      = CDF_TT2000_from_UTC_parts(YEAR_END, MONTH_END, DATE_END, HOUR_END, MINUTE_END, 0.0,  0.0,  0.0,  0.0);
    //**************************************************************************************************************
    


  CDF__TIME_TT2000 epoch_padValue = CDF_TT2000_from_UTC_parts(9999.0,12.0,31.0,23.0,59.0,59.0,999.0,999.0,999.0);
  double           F_sum1,eta_sum1,eta_dB1;
  double           F_sum2,eta_sum2,eta_dB2;
  double           F_sum3,eta_sum3,eta_dB3;
  double           X,Y,Z;
  ///////////////////////////////////
  double           F2[data_sum][64];   //    **************************
    
    
  char             date[50];

  //////////////////////////////////
  // ************************************
  sprintf(date,"%.0f/%.0f/%.0f",YEAR_START,MONTH_START,DATE_START);
  // ***************************************
  // CDFディレクトリパスの設定
  // ***************************************
  if((path = getenv("CDF_DATA_PATH")) == NULL){
      printf("Not find the File\n");
  }
  Y_ = M_ = D_ = 0;
    for(CDF__TIME_TT2000 epoch = epoch_start;epoch <= epoch_end;){
        CDF_TT2000_to_UTC_parts(epoch,&year,&month,&day,&hour,&minute,&second,TT2000NULL);
        YY  = year - 2000.0;
        MM  = month;
        DD  = day;
        hh  = hour;
        mm  = minute;
        ss  = second;
        if((Y_ != YY) || (M_ != MM) || (D_ != DD))
        {
            Y_ = YY; M_ = MM; D_ = DD;
            // ***************************************
            // CDFファイルを閉じる
            // ***************************************
            if(id1 != NULL){
                CDF__Close(&id1);
                id1 = NULL;
            }
            if(id2 != NULL){CDF__Close(&id2);id2 = NULL;}
            if(id3 != NULL){CDF__Close(&id3);id3 = NULL;}
            // ***************************************
            // CDFファイルパスの設定
            // ***************************************
            sprintf(dirname1 ,"%s/MRR1/MRR_L" CDF_DATA_LEVEL "/%04d%02d",path,2000 + YY,MM);
            sprintf(dirname2 ,"%s/MRR2/MRR_L" CDF_DATA_LEVEL "/%04d%02d",path,2000 + YY,MM);
            sprintf(dirname3 ,"%s/MRR3/MRR_L" CDF_DATA_LEVEL "/%04d%02d",path,2000 + YY,MM);
            sprintf(filename1,"%s/%s_l"  CDF_DATA_LEVEL "_%s_%04d%02d%02d_v" CDF_DATA_VERSION CDF_DATA_SUFFIX ".cdf",dirname1,"mrr1",DH1,2000 + YY,MM,DD);
            sprintf(filename2,"%s/%s_l"  CDF_DATA_LEVEL "_%s_%04d%02d%02d_v" CDF_DATA_VERSION CDF_DATA_SUFFIX ".cdf",dirname2,"mrr2",DH2,2000 + YY,MM,DD);
            sprintf(filename3,"%s/%s_l"  CDF_DATA_LEVEL "_%s_%04d%02d%02d_v" CDF_DATA_VERSION CDF_DATA_SUFFIX ".cdf",dirname3,"mrr3",DH3,2000 + YY,MM,DD);
            // ***************************************
            // MRR L1 CDFファイルを開く
            // ***************************************
            if(CDF__Open(&id1,filename1) == -1)
            {
                id1 = NULL;
                printf("File not found\n\n\n");
            }
            if(CDF__Open(&id2,filename2) == -1)id2 = NULL;
            if(CDF__Open(&id3,filename3) == -1)id3 = NULL;
        }
        if((id1 == NULL) || (id2 == NULL) || (id3 == NULL))goto EXIT;
        // ***************************************
        // MRR L1 データの取得
        // ***************************************
        recNum = (3600 * hh + 60 * mm + ss) / 10;
        CDF_L1_MRR__Get_Data(&mrr1,id1,recNum);
        CDF_L1_MRR__Get_Data(&mrr2,id2,recNum);
        CDF_L1_MRR__Get_Data(&mrr3,id3,recNum);
        // ***************************************
        // 1レコードについての出力処理
        // ***************************************
        /////////////////////////////
        for(i=0;i<64;i++){
            if (mrr2.F[FNUM][i]>1.0e-10) //printf("Zero ERROR\n");
                 F2[count][i] = 10.0*log10(mrr2.F[FNUM][i]);
            else F2[count][i] = 10.0*log10(1.0e-10);          // mrr2
        }
            //printf("Reading\n");
        count ++;
        // ***************************************
        // DT秒後のepochを得る
        // ***************************************
        EXIT:
        epoch = CDF_TT2000_from_UTC_parts(year,month,day,hour,minute,second + DT,TT2000END);
    }
    
    double dh,ha,dv;
    //高度分解能 : dh
    sscanf(DH2,"%lf",&dh);
    
    //測定高度 : ha
    ha=dh*(double)FNUM;
    printf("ha : %lf\n",ha);      //          **************
    //速度分解能 : dv
    dv = mrr2.dV;
    
    //main process
    int     Vnum;       //速度のナンバー
    double  g=0;        //
    double  g_stan;   //gを規格化
    double  f_stan;   //fを規格化
    
    double  Ft[data_sum][300];     //F(t) のデータを持つ配列
    
    double vel, tau;
    double v_max,   v_min;
    double tau_max, tau_min;
    v_max = V_max;
    v_min = V_min;
    tau_min = (H0-ha)/v_max + T0;
    tau_max = (H0-ha)/v_min + T0;
    
    
    printf(" tau max min = %lf, %lf \n",tau_max, tau_min);
    
// ***************** plot v_min  v_max ************************************
    gp = popen("gnuplot -persist","w");
    fprintf(gp, "pl '-' w l\n");
    for(i=(int)tau_min; i<(int)tau_max; i++){
        fprintf(gp, "%d %lf\n",i,fx(i, ha, H0, T0));
    }
    pclose(gp);
// ************************************************************************
    
    double h0, dh0;
    int H0_imax = 200;
    double h0_Height[H0_imax+1];
    dh0 = 10.0;         //d h0 [m]
    h0 = dh0 *  (int)( (H0-1000.0) / dh0);
    for (int m = 0; m < H0_imax; m++ ){
        h0_Height[m] = h0;
        for(i=0; i<(data_sum-(tau_max-tau_min)/10); i++){
            g=0.0;
            g_stan=0.0;
            f_stan=0.0;
        
            //Tau loop
        
            for(j=0; j<(tau_max-tau_min)/10; j++){
                tau = j* 10.0;

            //v loop
                for(k=0; k<64; k++){
                    vel = k * dv;
                    g +=F2[i+j][k] * Zc(tau, vel, ha, h0, dv);
                    f_stan += pow(Zc(tau, vel, ha, h0, dv), 2);
                    g_stan += pow(F2[i+j][k], 2);
                }
            }
            if ( (f_stan <=1.0e-10) || (g_stan <=1.0e-10) ) printf("ERROR\n");
            Ft[i][m] = g / (sqrt(g_stan) * sqrt(f_stan));
        }
        h0 += dh0;    //   +10[m]
    }
    //********************************************************************************
    //F_mrr2.dat output
    fp = fopen(DIR"F_mrr2"DATE".dat" , "w");
    h=HOUR_START;
    m=MINUTE_START;
    for(j=0;j<64;j++){
        for(i=0;i<data_sum;i++){
            //  x time
            //if(MINUTE_START-10<0)fprintf(fp, "%.0f:0%.0f ",HOUR_START,MINUTE_START);
            //else fprintf(fp, "%.0f:%.0f ",HOUR_START,MINUTE_START);
            fprintf(fp, "%.0f:%.0f:%.0f ",h,m,se);
            //  y dpsp
            /*
            if(i==0){
                fprintf(fp, "0 ");
            }else{*/
                fprintf(fp, "%f ",mrr2.V[j]);
           // }
            //  z power
            fprintf(fp, "%lf\n",F2[i][j]);
            //time*************************************
            se +=10;
            if( (se==60) & (m==59) ){
                h++;
                m=0;
                se=0;
            }else if(se==60){
                m++;
                se=0;
            }
            //******************************************
        }
        fprintf(fp, "\n");
        //reset
        h=HOUR_START;
        m=MINUTE_START;
        //
    }
    fclose(fp);
// **************************************************************************************
    
    //Ft(x).dat output

    fp = fopen(DIR"Ft_"DATE".dat" , "w");
    for(int n=0; n<H0_imax; n+=2){
        se=0;
        h=HOUR_START;
        m=MINUTE_START;
    //d[1]=2000+dh;  d[2]=240.4;
    for(i=0;i<(data_sum-(tau_max-tau_min)/10);i++){
        //  x time
        //if(MINUTE_START-10<0)
        //    fprintf(fp, "%.0f:%.0f:%.1f ",HOUR_START,MINUTE_START,se);
        //else
        //    fprintf(fp, "%.0f:%.0f:%.1f ",HOUR_START,MINUTE_START,se);
             fprintf(fp, "%.0f:%.0f:%.0f ",h,m,se);
        //  y dpsp
        fprintf(fp, "%lf %lf\n",h0_Height[n],Ft[i][n]);
        
        //time*************************************
        se +=10;
        if( (se>=60) & (m>=59) ){
            h++;
            m=0;
            se=0;
        }else if(se>=60){
            m++;
            se=0;
        }
            //******************************************
    
    }
        fprintf(fp, "\n");
        h=HOUR_START;
        m=MINUTE_START;
    }
    fclose(fp);
// *****************************************************************************************
    
    //********************************
    //gnuplot
    //********************************
    gp = popen("gnuplot -persist","w");
    
    fprintf(gp, "set terminal png size 1200,800  font 'Times New Roman, 18'\n");
    fprintf(gp, "unset key\n");

    fprintf(gp, "set xlabel 'Time' \n");
    //  set xadata
    fprintf(gp, "set xdata time\n");
    fprintf(gp, "set timefmt '%%H:%%M:%%S'\n");
    fprintf(gp, "set format x '%%H:%%M'\n");
    fprintf(gp, "set xtics 600\n");
    //  set graph
    if(m==0){
        fprintf(gp, "set xrange['%.0f:00':'%.0f:59']\n",h,HOUR_END-1);
    }else{
        fprintf(gp, "set xrange['%.0f:%.0f':'%.0f:%.0f']\n",h,m,HOUR_END,MINUTE_END-8-1);
    }
    
    fprintf(gp, "set origin 0, 0.03\n");
    fprintf(gp, "set xlabel font 'Times New Roman, 17'\n");
    fprintf(gp, "set ylabel font 'Times New Roman, 20'\n");
    fprintf(gp, "set xtics font 'Times New Roman, 14'\n");
    fprintf(gp, "set ytics font 'Times New Roman, 14'\n");
    //fprintf(gp, "set xlabel offset 0.5,-0.8\n");
    //fprintf(gp, "set ylabel offset 0,0\n");
    
    //2D graph : F(t)
    fprintf(gp, "set output '%sFt_%s.png' \n",DIR,DATE);
    fprintf(gp, "plot  '%sFt_%s.dat' using  1:3 w l linewidth 3\n",DIR,DATE);

    
    //3D graph
    
    fprintf(gp, "set ylabel 'Hight [m]' \n");
    fprintf(gp, "set yrange[%f:%f]\n",h0_Height[0],h0_Height[H0_imax-1]);
    if(m==0){
        fprintf(gp, "set xrange['%.0f:00':'%.0f:59']\n",h,HOUR_END-1);
    }else{
        fprintf(gp, "set xrange['%.0f:%.0f':'%.0f:%.0f']\n",h,m,HOUR_END,MINUTE_END-1);
    }
    //fprintf(gp, "set multiplot\n");
    fprintf(gp, "set pm3d map\n");
//    fprintf(gp, "set palette defined (0 'white', 1 'black')\n");
    
    //fprintf(gp, "set palette defined (0 'white', 0.3 'black', 1 'black')\n");
//    fprintf(gp, "set logscale cb\n");
    fprintf(gp, "set cbrange[0.14:0.2]\n");
    fprintf(gp, "set ztics font 'Times New Roman, 16'\n");
    
    
    fprintf(gp, "set output '%sFt_3d_%s.png' \n",DIR,DATE);
    //printf(gp, "set title 'F mrr2 %s (%d[m])' \n",date,mrr2.H[FNUM]);
    fprintf(gp, "splot  '%sFt_%s.dat' using  1:2:3\n",DIR,DATE);

 // **************************************************
    
    fprintf(gp, "set ylabel 'Velocity [m/s]' \n");
    fprintf(gp, "set yrange[0:%lf]\n",mrr2.V[63]);
    if(m==0){
        fprintf(gp, "set xrange['%.0f:00':'%.0f:59']\n",h,HOUR_END-1);
    }else{
        fprintf(gp, "set xrange['%.0f:%.0f':'%.0f:%.0f']\n",h,m,HOUR_END,MINUTE_END-1);
    }
    //fprintf(gp, "set multiplot\n");
    fprintf(gp, "set pm3d map\n");
    //fprintf(gp, "set palette defined (0 'white', 1 'black')\n");
    fprintf(gp, "set pm3d interpolate 5,5\n");
    //fprintf(gp, "set palette defined (0 'white', 0.3 'black', 1 'black')\n");
    fprintf(gp, "set cbrange[10:30.0]\n");
    fprintf(gp, "set ztics font 'Times New Roman, 16'\n");
    

    fprintf(gp, "set output '%sF_mrr2_%s.png' \n",DIR,DATE);
    //printf(gp, "set title 'F mrr2 %s (%d[m])' \n",date,mrr2.H[FNUM]);
    fprintf(gp, "splot  '%sF_mrr2%s.dat' using  1:2:3\n",DIR,DATE);
    
    
    pclose(gp);
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    //**********************************************************
    //**********************************************************

    /////////////////////////////////////////////////////////////////////////////////////////
    
    
    
  // ***************************************
  // CDFファイルを閉じる
  // ***************************************
    if(id1 != NULL)
    {
        CDF__Close(&id1);
        id1 = NULL;
    }
    if(id2 != NULL){CDF__Close(&id2);id2 = NULL;}
    if(id3 != NULL){CDF__Close(&id3);id3 = NULL;}

  return 0;
}


double  fx(double tau, double ha, double h0, double t0){
    double fi=0.0,ti=0.0;
    ti = tau - t0;
    if (ti != 0.0)
      fi = ( h0 - ha ) / ti;
    else printf(" NA overflow \n");
    return fi;
}

double Zc(double tau, double vel, double ha, double h0, double dv)
{
    double t0;
    t0 = T0;
    double vel_t;
    vel_t = fx(tau, ha, h0, t0);
    if (  ((vel_t-dv/2.0) <= vel )  && ((vel_t+dv/2.0)> vel ))
        return 1.0;
    else
        return 0.0;
}

