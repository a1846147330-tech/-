// V5：统一16变量凸QP，利用辅助变量块结构将Newton分解降到8维。
// 不依赖第三方数值内核；C++17标准库与Windows API，全部计算在独立进程内完成。
// 阅读顺序：make_qp -> solve_qp -> optimize -> calculate -> main。
// 组限制和 L1 换手预算来自题目 algo.py 的 get_constrain，并非额外添加。
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <array>
#include <cstdint>
#include <cstring>



#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>

struct InputData {
    std::vector<int> dates;  // Excel 1900-system serial date, whole days
    std::vector<std::array<double, 8>> values;
};


inline InputData read_xls(const std::string& path);
#include <algorithm>
#include <charconv>
#include <chrono>




#include <limits>
#include <numeric>

// 【参数】N=8 个策略；WINDOW=60 行历史；KEEP=15 行持有期。
constexpr int N=8, M=18, WINDOW=60, KEEP=15;
// Vec 类似 Python 中长度为 8 的列表；Series 类似二维列表。
using Vec=std::array<double,N>;
using Matrix=std::array<Vec,N>;
using Series=std::vector<Vec>;
const Vec INITIAL={8.9/100.,10.25/100.,6.35/100.,6.45/100.,25.72/100.,16.55/100.,12.72/100.,13.06/100.};
struct Moments { Vec mean{},numerator{}; Matrix cov{}; };
struct SolveResult { Vec weight; int status,iterations; double objective,violation; double residual=0; };
struct RunResult {
    std::string id;
    Series weights;
    std::vector<double> returns;
    std::vector<SolveResult> solves;
};

inline double dot(const Vec&a,const Vec&b){
    double result=0;for(int j=0;j<N;++j)result+=a[j]*b[j];return result;
}
inline double sum8(const Vec&a){return ((a[0]+a[1])+(a[2]+a[3]))+((a[4]+a[5])+(a[6]+a[7]));}
template<class Container> inline double pairwise_sum(const Container&a){
    // NumPy's small contiguous reduction uses eight partial accumulators.
    if(a.size()<8){double result=-0.;for(double x:a)result+=x;return result;}
    Vec sums;std::copy_n(a.begin(),8,sums.begin());size_t i=8;
    for(;i+7<a.size();i+=8)for(int j=0;j<8;++j)sums[j]+=a[i+j];
    double result=sum8(sums);for(;i<a.size();++i)result+=a[i];return result;
}
inline Vec multiply(const Matrix&a,const Vec&b){
    Vec out{};for(int j=0;j<N;++j)for(int k=0;k<N;++k)out[j]+=a[j][k]*b[k];return out;
}
// 【1. 差分】对应 Python：df.diff(1).fillna(0)。
inline Series differences(const Series& values){
    Series r(values.size());
    for(size_t i=1;i<r.size();++i)for(int j=0;j<N;++j)r[i][j]=values[i][j]-values[i-1][j];
    return r;
}
// 【2. 去极值】对应 norm_preprocess；处理副本，不覆盖原始收益。
inline Series preprocess(Series p){
    // Deliberately use earlier *clipped* rows, as algo.py does.
    for(size_t i=WINDOW;i<p.size();++i)for(int j=0;j<N;++j){
        // 固定 60 点工作区放在栈上；保留原来的分组求和顺序。
        std::array<double,WINDOW> values;for(int k=0;k<WINDOW;++k)values[k]=p[i-WINDOW+k][j];
        double avg=pairwise_sum(values)/WINDOW;
        for(double&v:values)v=(avg-v)*(avg-v);
        double std=std::sqrt(pairwise_sum(values)/(WINDOW-1));
        p[i][j]=std::clamp(p[i][j],avg-3*std,avg+3*std);
    }
    return p;
}
// 【3. 统计】先取 60 行，再滚动求 15 行之和，共 46 个样本。
// mean 类似长度 8 的 numpy 数组，cov 类似 8×8 数组。
inline Moments moments(const Series&p,size_t end,bool ewm){
    // 60 rows -> 46 overlapping 15-row sums. No observations from end or later.
    std::array<Vec,WINDOW-KEEP+1> sums{};
    for(int j=0;j<N;++j){
        double total=0,add_comp=0,remove_comp=0;
        for(int k=0;k<WINDOW;++k){
            if(k>=KEEP){double y=-p[end-WINDOW+k-KEEP][j]-remove_comp;
                double t=total+y;remove_comp=(t-total)-y;total=t;}
            double y=p[end-WINDOW+k][j]-add_comp,t=total+y;add_comp=(t-total)-y;total=t;
            if(k>=KEEP-1)sums[k-KEEP+1][j]=total;
        }
    }
    Moments m;
    if(!ewm){
        for(int j=0;j<N;++j){
            std::array<double,WINDOW-KEEP+1> column;for(size_t i=0;i<sums.size();++i)column[i]=sums[i][j];
            m.mean[j]=pairwise_sum(column)/sums.size();
            for(double&v:column)v=(m.mean[j]-v)*(m.mean[j]-v);
            m.numerator[j]=std::sqrt(pairwise_sum(column)/(sums.size()-1));
        }
        int n=N,k=int(sums.size());std::array<double,N*(WINDOW-KEEP+1)> centered;
        for(int j=0;j<n;++j)for(int i=0;i<k;++i)centered[j*k+i]=sums[i][j]-m.mean[j];
        for(int j=0;j<n;++j)for(int i=0;i<=j;++i){
            double value=0;for(int r=0;r<k;++r)value+=centered[j*k+r]*centered[i*k+r];
            m.cov[j][i]=m.cov[i][j]=value/(k-1);
        }
    }else{
        // pandas adjust=False, bias=False: first observation has weight 1.
        const double alpha=2./(WINDOW+1), beta=1-alpha;
        m.mean=sums[0];double sumw2=1;
        for(size_t i=1;i<sums.size();++i){
            Vec old=m.mean;
            for(int j=0;j<N;++j)if(old[j]!=sums[i][j])m.mean[j]=beta*old[j]+alpha*sums[i][j];
            for(int j=0;j<N;++j)for(int k=0;k<N;++k)
                m.cov[j][k]=beta*(m.cov[j][k]+(old[j]-m.mean[j])*(old[k]-m.mean[k]))
                    +alpha*((sums[i][j]-m.mean[j])*(sums[i][k]-m.mean[k]));
            sumw2=beta*beta*sumw2+alpha*alpha;
        }
        for(auto&r:m.cov)for(double&v:r)v*=1./(1-sumw2);
        // Compatibility: source uses variance, not sqrt(variance), in this branch.
        for(int j=0;j<N;++j)m.numerator[j]=m.cov[j][j];
    }
    return m;
}
// 【4. 目标函数】risk<0 为最大分散；risk=0/20 为均值方差。
inline double objective(const Vec&w,const Moments&m,double risk){
    double q=dot(multiply(m.cov,w),w);
    if(risk<0){
        if(q<=0)throw std::runtime_error("Nonpositive portfolio variance");
        return -dot(w,m.numerator)/std::sqrt(q);
    }
    Vec prod;for(int j=0;j<N;++j)prod[j]=w[j]*m.mean[j];
    return -sum8(prod)+.5*risk*q;
}
// 【5. 约束】第 0 项须等于 0，其余项须 >=0；对应 get_constrain。
inline std::array<double,M> constraints(const Vec&w,const Vec&prev){
    std::array<double,M> c{};
    c[0]=sum8(w)-1;
    for(int j=0;j<N;++j)c[j+1]=w[j];
    Vec changes;for(int j=0;j<N;++j)changes[j]=std::abs(w[j]-prev[j]);c[9]=.1-sum8(changes);
    const double groups[]={w[0]+w[1]+w[2]+w[3],w[4],w[5],w[6]+w[7]};
    for(int k=0;k<4;++k){c[10+2*k]=groups[k]-.15;c[11+2*k]=.35-groups[k];}
    return c;
}
inline double constraint_violation(const Vec&w,const Vec&prev){
    auto c=constraints(w,prev);double v=std::abs(c[0]);
    for(int i=1;i<M;++i)v=std::max(v,-c[i]);return v;
}

// 【核心创新】最大分散比率的正齐次变换，使三类目标共用一个小规模凸 QP。
// 变量 x=(u,t)：8 个组合变量 + 8 个换手辅助变量。
// 普通目标 u=w，sum(w)=1；分散目标 u=w/(a'w)，a'u=1，w=u/sum(u)。
// 后者的所有组限制和换手预算随 sum(u) 一起缩放，可行集完全等价。
constexpr int NV=16, NC=33;
using V16=std::array<double,NV>;
struct QP {
    std::array<V16,NV> H{};
    V16 q{},e{},initial{};
    std::array<V16,NC> A{};
    std::array<double,NC> b{};
    std::array<std::array<int,NV>,NC> nz{};
    std::array<int,NC> nnz{};
};
inline QP make_qp(const Vec&prev,const Moments&m,double risk){
    QP qp;bool div=risk<0;
    Vec a=m.numerator;
    double amax=*std::max_element(a.begin(),a.end());
    if(div && !(amax>0))throw std::runtime_error("No positive diversification numerator");
    if(div)for(double&v:a)v/=amax; // 归一化分子只改变比率尺度，不改变最优权重。
    double scale=div?dot(a,prev):1.;
    for(int j=0;j<N;++j){qp.e[j]=div?a[j]:1;qp.initial[j]=prev[j]/scale;qp.initial[N+j]=.05/N/scale;}
    double objective_scale=0;
    for(int j=0;j<N;++j){
        qp.q[j]=div?0:-m.mean[j];
        for(int k=0;k<N;++k)qp.H[j][k]=(div?2:risk)*m.cov[j][k];
        objective_scale=std::max(objective_scale,std::abs(qp.q[j]));
        for(int k=0;k<N;++k)objective_scale=std::max(objective_scale,std::abs(qp.H[j][k]));
    }
    if(objective_scale==0)objective_scale=1;
    for(int j=0;j<NV;++j){qp.q[j]/=objective_scale;for(double&v:qp.H[j])v/=objective_scale;}
    int row=0;
    for(int j=0;j<N;++j)qp.A[row++][j]=-1; // 非负权重。
    for(int g=0;g<4;++g){
        for(int j=0;j<N;++j){
            bool member=(g==0?j<4:g==1?j==4:g==2?j==5:j>=6);
            qp.A[row][j]=-(member?1.:0)+(div?.15:0);
            qp.A[row+1][j]=(member?1.:0)-(div?.35:0);
        }
        qp.b[row]=div?0:-.15;qp.b[row+1]=div?0:.35;row+=2;
    }
    for(int j=0;j<N;++j){
        for(int k=0;k<N;++k){qp.A[row][k]=(j==k?1.:0)-(div?prev[j]:0);qp.A[row+1][k]=-(j==k?1.:0)+(div?prev[j]:0);}
        qp.A[row][N+j]=qp.A[row+1][N+j]=-1;
        qp.b[row]=div?0:prev[j];qp.b[row+1]=div?0:-prev[j];row+=2;
    }
    for(int j=0;j<N;++j){qp.A[row][N+j]=1;qp.A[row][j]=div?-.1:0;}
    qp.b[row]=div?0:.1;
    // 记录每行非零列：普通换手约束只有两项，避免按 16×16 稠密矩阵累加零。
    for(int i=0;i<NC;++i)for(int j=0;j<NV;++j)if(qp.A[i][j]!=0)qp.nz[i][qp.nnz[i]++]=j;
    return qp;
}

struct QPSolution {V16 x;int iterations;double residual;};
// 小矩阵使用固定长度栈数组；同一次迭代的预测步与校正步共享 Cholesky 分解。
inline QPSolution solve_qp(const QP&p){
    V16 x=p.initial,rd{},dx{},affx{},ed{};
    std::array<double,NC> s{},z{},rp{},ds{},dz{},affs{},affz{},rc{};
    double nu=0,dnu=0;
    auto rowdot=[&](int i,const V16&v){double sum=0;for(int h=0;h<p.nnz[i];++h){int j=p.nz[i][h];sum+=p.A[i][j]*v[j];}return sum;};
    for(int i=0;i<NC;++i){s[i]=std::max(1.,p.b[i]-rowdot(i,x));z[i]=1;}
    double residual=1e100;
    for(int iter=0;iter<80;++iter){
        double re=-1,mu=0;
        for(int j=0;j<NV;++j){re+=p.e[j]*x[j];rd[j]=p.q[j]+p.e[j]*nu;
            if(j<N)for(int k=0;k<N;++k)rd[j]+=p.H[j][k]*x[k];}
        for(int i=0;i<NC;++i)for(int h=0;h<p.nnz[i];++h){int j=p.nz[i][h];rd[j]+=p.A[i][j]*z[i];}
        residual=std::abs(re);
        for(int i=0;i<NC;++i){rp[i]=s[i]-p.b[i]+rowdot(i,x);
            residual=std::max(residual,std::abs(rp[i]));mu+=s[i]*z[i];}
        mu/=NC;residual=std::max(residual,mu);
        for(double v:rd)residual=std::max(residual,std::abs(v));
        if(residual<1e-9)return {x,iter,residual};
        std::array<V16,NV> L=p.H;
        for(int i=0;i<NC;++i){double ratio=z[i]/s[i];
            for(int h=0;h<p.nnz[i];++h){int j=p.nz[i][h];double aj=ratio*p.A[i][j];
                for(int l=0;l<=h;++l){int k=p.nz[i][l];L[j][k]+=aj*p.A[i][k];}}}
        // 微小正则项仅作用于 Newton 线性系统，收敛判断仍针对未正则化的原问题。
        for(int j=0;j<NV;++j)L[j][j]+=1e-10;
        // 辅助换手变量块 T=diag(d)+rho*11'；Sherman-Morrison 消元后仅分解8×8。
        // 若相减损失相对精度或 Schur 补无法稳定分解，回退原16×16系统。
        // 原K矩阵L保留至快路径成功，原问题的1e-9残差判据始终不变。
        Matrix B{},S{};Vec dinv{},bd{};
        double rho=L[N+1][N],dsum=0,rank=0;
        bool reduced=std::isfinite(rho) && rho>=0;
        constexpr double guard=64*std::numeric_limits<double>::epsilon();
        for(int j=0;j<N && reduced;++j){
            double diagonal=L[N+j][N+j],d=diagonal-rho;
            reduced=std::isfinite(d) && d>guard*std::max(std::abs(diagonal),std::abs(rho));
            if(reduced){dinv[j]=1/d;dsum+=dinv[j];}
        }
        if(reduced){double denominator=1+rho*dsum;
            reduced=std::isfinite(dsum) && std::isfinite(denominator) && denominator>0;
            if(reduced){rank=rho/denominator;reduced=std::isfinite(rank);}}
        if(reduced){
            for(int i=0;i<N;++i)for(int j=0;j<N;++j){B[i][j]=L[N+j][i];bd[i]+=B[i][j]*dinv[j];}
            for(int i=0;i<N;++i)for(int j=0;j<=i;++j){
                double v=L[i][j];for(int k=0;k<N;++k)v-=B[i][k]*B[j][k]*dinv[k];
                S[i][j]=v+rank*bd[i]*bd[j];
            }
            for(int j=0;j<N && reduced;++j)for(int k=0;k<=j && reduced;++k){
                double value=S[j][k];for(int h=0;h<k;++h)value-=S[j][h]*S[k][h];
                if(j==k){
                    reduced=std::isfinite(value) && value>guard*std::abs(S[j][j]);
                    if(reduced)S[j][k]=std::sqrt(value);
                }else {S[j][k]=value/S[k][k];reduced=std::isfinite(S[j][k]);}
            }
        }
        if(!reduced){
            for(int j=0;j<NV;++j)for(int k=0;k<=j;++k){double value=L[j][k];for(int h=0;h<k;++h)value-=L[j][h]*L[k][h];
                if(j==k){if(!(value>0) || !std::isfinite(value))throw std::runtime_error("Newton matrix is not positive definite");L[j][k]=std::sqrt(value);}
                else L[j][k]=value/L[k][k];}
        }
        auto tinv=[&](Vec v){double d=0;for(int j=0;j<N;++j){v[j]*=dinv[j];d+=v[j];}
            d*=rank;for(int j=0;j<N;++j)v[j]-=dinv[j]*d;return v;};
        auto solve=[&](V16 v){
            if(!reduced){
                for(int j=0;j<NV;++j){for(int k=0;k<j;++k)v[j]-=L[j][k]*v[k];v[j]/=L[j][j];}
                for(int j=NV-1;j>=0;--j){for(int k=j+1;k<NV;++k)v[j]-=L[k][j]*v[k];v[j]/=L[j][j];}return v;
            }
            Vec t;for(int j=0;j<N;++j)t[j]=v[N+j];t=tinv(t);
            for(int j=0;j<N;++j)for(int k=0;k<N;++k)v[j]-=B[j][k]*t[k];
            for(int j=0;j<N;++j){for(int k=0;k<j;++k)v[j]-=S[j][k]*v[k];v[j]/=S[j][j];}
            for(int j=N-1;j>=0;--j){for(int k=j+1;k<N;++k)v[j]-=S[k][j]*v[k];v[j]/=S[j][j];}
            for(int j=0;j<N;++j){t[j]=v[N+j];for(int k=0;k<N;++k)t[j]-=B[k][j]*v[k];}
            t=tinv(t);for(int j=0;j<N;++j)v[N+j]=t[j];return v;
        };

        ed=solve(p.e);double denom=0;for(int j=0;j<NV;++j)denom+=p.e[j]*ed[j];
        auto direction=[&](){
            V16 rhs;for(int j=0;j<NV;++j)rhs[j]=-rd[j];
            for(int i=0;i<NC;++i){double value=(rc[i]-z[i]*rp[i])/s[i];for(int h=0;h<p.nnz[i];++h){int j=p.nz[i][h];rhs[j]+=p.A[i][j]*value;}}
            dx=solve(rhs);double eq=0;for(int j=0;j<NV;++j)eq+=p.e[j]*dx[j];dnu=(eq+re)/denom;
            for(int j=0;j<NV;++j)dx[j]-=ed[j]*dnu;
            for(int i=0;i<NC;++i){ds[i]=-rp[i]-rowdot(i,dx);dz[i]=(-rc[i]-z[i]*ds[i])/s[i];}
        };
        auto steps=[&](){std::array<double,2> step={1,1};for(int i=0;i<NC;++i){if(ds[i]<0)step[0]=std::min(step[0],-s[i]/ds[i]);if(dz[i]<0)step[1]=std::min(step[1],-z[i]/dz[i]);}return step;};
        for(int i=0;i<NC;++i)rc[i]=s[i]*z[i];direction();auto alpha=steps();
        double mua=0;for(int i=0;i<NC;++i){mua+=(s[i]+alpha[0]*ds[i])*(z[i]+alpha[1]*dz[i]);affs[i]=ds[i];affz[i]=dz[i];}
        double sigma=std::pow(std::clamp(mua/(NC*mu),0.,1.),3);
        for(int i=0;i<NC;++i)rc[i]=s[i]*z[i]+affs[i]*affz[i]-sigma*mu;
        direction();alpha=steps();alpha[0]=std::min(1.,.995*alpha[0]);alpha[1]=std::min(1.,.995*alpha[1]);
        for(int j=0;j<NV;++j)x[j]+=alpha[0]*dx[j];nu+=alpha[1]*dnu;
        for(int i=0;i<NC;++i){s[i]+=alpha[0]*ds[i];z[i]+=alpha[1]*dz[i];}
    }
    throw std::runtime_error("Interior point method did not converge; residual="+std::to_string(residual));
}
inline SolveResult optimize(const Vec&prev,const Moments&m,double risk){
    auto sol=solve_qp(make_qp(prev,m,risk));Vec w;
    double scale=1;if(risk<0){scale=0;for(int j=0;j<N;++j)scale+=sol.x[j];}
    for(int j=0;j<N;++j){w[j]=sol.x[j]/scale;if(!std::isfinite(w[j]))throw std::runtime_error("Nonfinite weight");}
    double violation=constraint_violation(w,prev);
    if(violation>1e-7)throw std::runtime_error("Weight constraints failed");
    return {w,0,sol.iterations,objective(w,m,risk),violation,sol.residual};
}

// 【7. 回测主循环】最值得先看：对应 get_weight + run。
// 每种配置调仓 12 次；每次用 46 个样本求一组权重，不是求 46 次。
inline std::vector<RunResult> calculate(const InputData&data){
    auto r=differences(data.values),p=preprocess(r);
    // 三种目标使用相同历史统计；每个调仓日只计算两种估计，各一次。
    std::vector<std::array<Moments,2>> history;
    for(size_t i=WINDOW;i<r.size();i+=KEEP)
        history.push_back({moments(p,i,false),moments(p,i,true)});
    std::vector<RunResult> results;
    const double risks[]={-1,0,20};const char* names[]={"div","ret","mvo"};
    for(int mode=0;mode<3;++mode)for(bool ewm:{false,true}){
        RunResult result;result.id=std::string(names[mode])+(ewm?"_ewm":"_equal");
        result.weights.reserve(r.size());
        result.returns.reserve(r.size());
        result.solves.reserve(history.size());
        Vec w=INITIAL;
        for(size_t i=0;i<r.size();++i){
            if(i>=WINDOW && (i-WINDOW)%KEEP==0){
                auto solution=optimize(w,history[(i-WINDOW)/KEEP][ewm?1:0],risks[mode]);
                result.solves.push_back(solution);w=solution.weight;
            }
            result.weights.push_back(w);result.returns.push_back(dot(w,r[i]));
        }
        results.push_back(result);
    }
    RunResult base;base.id="baseline";
    for(const auto&v:r){base.weights.push_back(INITIAL);base.returns.push_back(dot(INITIAL,v));}
    results.push_back(base);return results;
}
inline std::string date_string(int serial){
    // Excel 日序号直接转公历；全整数运算，避免逐日调用时区和区域格式化。
    int z=serial-25569+719468;
    int era=(z>=0?z:z-146096)/146097;
    unsigned doe=unsigned(z-era*146097);
    unsigned yoe=(doe-doe/1460+doe/36524-doe/146096)/365;
    int y=int(yoe)+era*400;
    unsigned doy=doe-(365*yoe+yoe/4-yoe/100),mp=(5*doy+2)/153;
    unsigned d=doy-(153*mp+2)/5+1,m=mp<10?mp+3:mp-9;
    y+=(m<=2);
    char text[8]={char('0'+y/1000),char('0'+y/100%10),char('0'+y/10%10),char('0'+y%10),char('0'+m/10),char('0'+m%10),char('0'+d/10),char('0'+d%10)};
    return std::string(text,8);
}
struct Metrics{std::array<double,6> numbers;std::string peak,trough;};
// 【8. 评价指标】对应 EvaluatePNL。
inline Metrics evaluate(const RunResult&r,const InputData&data){
    std::vector<double> pnl(r.returns.size()),daily;double sum=0;
    for(size_t i=0;i<pnl.size();++i){sum+=r.returns[i];pnl[i]=sum+1;
        if(i && pnl[i]!=pnl[i-1]) daily.push_back(pnl[i]-pnl[i-1]);}
    if(daily.size()<2)throw std::runtime_error("Insufficient nonzero returns for metrics");
    double mean=std::accumulate(daily.begin(),daily.end(),0.)/daily.size(),var=0;
    size_t wins=0;for(double v:daily){var+=(v-mean)*(v-mean);wins+=(v>0);}
    double vol=std::sqrt(var/(daily.size()-1))*std::sqrt(250.),annual=mean*250;
    double dd=-std::numeric_limits<double>::infinity();size_t peak=0,left=0,right=1;
    for(size_t j=1;j<pnl.size();++j){
        double d=pnl[peak]-pnl[j];
        if(d>dd || (d==dd && (peak<left || (peak==left && j<right)))){dd=d;left=peak;right=j;}
        if(pnl[j]>pnl[peak])peak=j;
    }
    double years=(data.dates.back()-data.dates.front())/365.;
    return {{annual,vol,annual/vol,(pnl.back()-pnl.front())/(years*dd),double(wins)/daily.size(),dd},
        date_string(data.dates[left]),date_string(data.dates[right])};
}
// 【导出】单个 JSON 包含指标、全部日收益、权重及 72 次求解证据。
inline void number(std::string&out,double value){char b[64];auto r=std::to_chars(b,b+64,value,std::chars_format::general,17);if(r.ec!=std::errc{} || !std::isfinite(value))throw std::runtime_error("Invalid output number");out.append(b,r.ptr);}
inline void dump(const std::string&path,const InputData&data,const std::vector<RunResult>&runs){
    std::string out="{\"method\":\"convex_qp_ipm\",\"dates\":[";out.reserve(400000);
    for(size_t i=0;i<data.dates.size();++i){if(i)out+=',';out+='"'+date_string(data.dates[i])+'"';}out+="],\"models\":[";
    for(size_t k=0;k<runs.size();++k){auto&r=runs[k];auto m=evaluate(r,data);if(k)out+=',';
        out+="{\"id\":\""+r.id+"\",\"metrics\":[";for(int j=0;j<6;++j){if(j)out+=',';number(out,m.numbers[j]);}out+=",\""+m.peak+"\",\""+m.trough+"\"],\"returns\":[";
        for(size_t i=0;i<r.returns.size();++i){if(i)out+=',';number(out,r.returns[i]);}out+="],\"weights\":[";
        // 同一持有期权重不变：只复用本次已计算权重的文本，仍逐日完整写出。
        // 比较实际向量，既不硬编码调仓日，也不缓存跨请求的数值答案。
        std::string weight_text;
        for(size_t i=0;i<r.weights.size();++i){
            if(i)out+=',';
            if(i==0 || r.weights[i]!=r.weights[i-1]){
                weight_text.clear();weight_text+='[';
                for(int j=0;j<N;++j){if(j)weight_text+=',';number(weight_text,r.weights[i][j]);}
                weight_text+=']';
            }
            out+=weight_text;
        }out+="],\"solver\":[";
        for(size_t i=0;i<r.solves.size();++i){if(i)out+=',';auto&s=r.solves[i];out+="{\"iterations\":"+std::to_string(s.iterations)+",\"objective\":";number(out,s.objective);out+=",\"violation\":";number(out,s.violation);out+=",\"kkt_residual\":";number(out,s.residual);out+='}';}out+="]}";
    }out+="]}";
    // 直接写出缓冲并关闭句柄，不触发 iostream 的区域对象初始化。
    HANDLE f=CreateFileA(path.c_str(),GENERIC_WRITE,FILE_SHARE_READ,nullptr,CREATE_ALWAYS,FILE_ATTRIBUTE_NORMAL,nullptr);
    if(f==INVALID_HANDLE_VALUE)throw std::runtime_error("Cannot open JSON");
    DWORD written=0;BOOL ok=WriteFile(f,out.data(),DWORD(out.size()),&written,nullptr);BOOL closed=CloseHandle(f);
    if(!ok || written!=out.size() || !closed)throw std::runtime_error("Cannot write complete JSON");
}
inline void write_console(DWORD kind,const std::string&text){
    HANDLE out=GetStdHandle(kind);if(out==nullptr || out==INVALID_HANDLE_VALUE)return;
    DWORD written=0;WriteFile(out,text.data(),DWORD(text.size()),&written,nullptr);
}
// 独立结构检查：扰动未来净值不应影响此前权重；回撤与全枚举逐项比较。
inline void self_test(const InputData&data){
    auto base=calculate(data);auto changed=data;
    size_t split=data.values.size()*3/4;
    for(size_t i=split;i<changed.values.size();++i)for(int j=0;j<N;++j)changed.values[i][j]+=.003*(j+1);
    auto future=calculate(changed);
    for(size_t k=0;k<base.size();++k){
        for(size_t i=0;i<=split;++i)if(base[k].weights[i]!=future[k].weights[i])throw std::runtime_error("Lookahead test failed");
        auto metric=evaluate(base[k],data);std::vector<double> pnl;double level=1;
        for(double r:base[k].returns){level+=r;pnl.push_back(level);}
        double best=-1e100;size_t left=0,right=1;
        for(size_t i=0;i<pnl.size();++i)for(size_t j=i+1;j<pnl.size();++j)if(pnl[i]-pnl[j]>best){best=pnl[i]-pnl[j];left=i;right=j;}
        if(std::abs(best-metric.numbers[5])>1e-12 || date_string(data.dates[left])!=metric.peak || date_string(data.dates[right])!=metric.trough)throw std::runtime_error("Drawdown test failed");
    }
    write_console(STD_OUTPUT_HANDLE,"Self-test passed: future-data isolation and exhaustive drawdown comparison\n");
}
// 标准输出使用一次 WriteFile；无 iostream 静态初始化。
int main(int argc,char**argv){
    try{
        if(argc>2 && std::string(argv[1])=="--self-test"){self_test(read_xls(argv[2]));return 0;}
        auto start=std::chrono::steady_clock::now();
        char module[MAX_PATH];DWORD length=GetModuleFileNameA(nullptr,module,MAX_PATH);
        if(!length || length>=MAX_PATH)throw std::runtime_error("Executable path is too long");
        std::string exe(module,length);exe.resize(exe.find_last_of("\\/"));
        std::string input=argc>1?argv[1]:exe.substr(0,exe.find_last_of("\\/"))+"/data.xls";
        std::string output=argc>2?argv[2]:exe+"/results.json";
        int repeat=argc>3?std::stoi(argv[3]):1;if(repeat<1 || repeat>10000)throw std::runtime_error("repeat must be 1..10000");
        auto data=read_xls(input);auto loaded=std::chrono::steady_clock::now();
        std::vector<RunResult> runs;std::vector<double> samples;
        for(int i=0;i<repeat;++i){auto t=std::chrono::steady_clock::now();runs=calculate(data);for(auto&r:runs)(void)evaluate(r,data);samples.push_back(std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count());}
        auto computed=std::chrono::steady_clock::now();dump(output,data,runs);auto written=std::chrono::steady_clock::now();
        std::string text;text.reserve(2000);
        for(auto&r:runs){auto m=evaluate(r,data);text+=r.id;for(double v:m.numbers){text+=' ';number(text,v);}text+=' '+m.peak+' '+m.trough+'\n';}
        std::sort(samples.begin(),samples.end());auto ms=[](auto a,auto b){return std::chrono::duration<double,std::milli>(b-a).count();};
        text+="TIMING load_ms=";number(text,ms(start,loaded));text+=" compute_min_ms=";number(text,samples.front());text+=" compute_median_ms=";number(text,samples[samples.size()/2]);text+=" compute_max_ms=";number(text,samples.back());text+=" export_ms=";number(text,ms(computed,written));text+=" total_ms=";number(text,ms(start,written));text+=" repeats="+std::to_string(repeat)+'\n';
        write_console(STD_OUTPUT_HANDLE,text);
        return 0;
    }catch(const std::exception&e){write_console(STD_ERROR_HANDLE,std::string("Error: ")+e.what()+'\n');return 1;}
}

// 底层 Excel 读取：本题限定格式，可暂时跳过。
inline InputData read_xls(const std::string& path) {
    HANDLE f=CreateFileA(path.c_str(),GENERIC_READ,FILE_SHARE_READ,nullptr,OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,nullptr);
    if(f==INVALID_HANDLE_VALUE)throw std::runtime_error("Cannot open XLS: "+path);
    LARGE_INTEGER length{};BOOL sized=GetFileSizeEx(f,&length);
    if(!sized || length.QuadPart<512 || length.QuadPart>64*1024*1024){CloseHandle(f);throw std::runtime_error("Invalid XLS size");}
    std::vector<unsigned char> file(size_t(length.QuadPart));
    DWORD read=0;BOOL ok=ReadFile(f,file.data(),DWORD(file.size()),&read,nullptr);BOOL closed=CloseHandle(f);
    if(!ok || read!=file.size() || !closed)throw std::runtime_error("Cannot read complete XLS");
    auto u16 = [](const std::vector<unsigned char>& b, size_t p) -> unsigned {
        if (p + 2 > b.size()) throw std::runtime_error("Truncated XLS record");
        return b[p] | unsigned(b[p+1]) << 8;
    };
    auto u32 = [&](const std::vector<unsigned char>& b, size_t p) -> uint32_t {
        return u16(b,p) | uint32_t(u16(b,p+2)) << 16;
    };
    const unsigned char signature[] = {0xd0,0xcf,0x11,0xe0,0xa1,0xb1,0x1a,0xe1};
    if (file.size()<512 || std::memcmp(file.data(),signature,8))
        throw std::runtime_error("Input must be an OLE BIFF8 .xls file");
    if (u16(file,30)!=9 || u16(file,28)!=0xfffe || u32(file,72)!=0)
        throw std::runtime_error("Only little-endian 512-byte CFB sectors / header DIFAT supported");
    auto sector = [&](uint32_t id) -> size_t {
        size_t p=(size_t(id)+1)*512;
        if(p+512>file.size()) throw std::runtime_error("Invalid CFB sector");
        return p;
    };
    std::vector<uint32_t> fat;
    unsigned nf=u32(file,44);
    if (nf>109) throw std::runtime_error("Extended DIFAT unsupported");
    for(unsigned i=0;i<nf;++i){
        size_t p=sector(u32(file,76+4*i));
        for(unsigned j=0;j<128;++j) fat.push_back(u32(file,p+4*j));
    }
    auto stream = [&](uint32_t id) {
        std::vector<unsigned char> b;
        std::vector<unsigned char> visited(fat.size());
        while(id!=0xfffffffe){
            if(id>=fat.size() || visited[id]++)
                throw std::runtime_error("Invalid/cyclic CFB chain");
            size_t p=sector(id);
            b.insert(b.end(),file.begin()+p,file.begin()+p+512);
            id=fat[id];
        }
        return b;
    };
    auto directory=stream(u32(file,48));
    std::vector<unsigned char> book;
    for(size_t p=0;p+128<=directory.size();p+=128){
        unsigned len=u16(directory,p+64);
        if(len<2 || len>64 || directory[p+66]!=2) continue;
        std::string name;
        for(unsigned j=0;j+2<len;j+=2) name+=char(u16(directory,p+j));
        if(name=="Workbook" || name=="Book"){
            uint32_t bytes=u32(directory,p+120);
            if(bytes<4096 || u32(directory,p+124)!=0)
                throw std::runtime_error("Mini/huge Workbook streams unsupported");
            book=stream(u32(directory,p+116));
            if(bytes>book.size()) throw std::runtime_error("Truncated Workbook stream");
            book.resize(bytes);
        }
    }
    if(book.empty()) throw std::runtime_error("Workbook stream absent");
    size_t sheet=0;
    for(size_t p=0;p+4<=book.size();){
        unsigned id=u16(book,p),len=u16(book,p+2);size_t d=p+4;
        if(d+len>book.size()) throw std::runtime_error("Truncated BIFF record");
        if(id==0x85 && !sheet){
            if(len<8 || book[d+5]!=0) throw std::runtime_error("First sheet is not a worksheet");
            sheet=u32(book,d);
        }
        if(id==0x22 && u16(book,d)!=0) throw std::runtime_error("1904 date system unsupported");
        if(id==0x2f) throw std::runtime_error("Encrypted XLS unsupported");
        if(id==0x0a) break;
        p=d+len;
    }
    if(!sheet) throw std::runtime_error("No sheet found");
    auto rk = [](uint32_t raw){
        double value;
        if(raw&2) value=double(int32_t(raw)>>2);
        else {uint64_t bits=uint64_t(raw&~3u)<<32;std::memcpy(&value,&bits,8);}
        return raw&1 ? value/100.0 : value;
    };
    std::vector<std::array<double,9>> rows(1);
    std::vector<unsigned> masks(1);
    auto put = [&](unsigned r,unsigned c,double value){
        if(r==0) return;
        if(c>8) throw std::runtime_error("Unexpected numeric column");
        if(r>=rows.size()){rows.resize(r+1);masks.resize(r+1);} 
        if(!std::isfinite(value) || (masks[r]&(1u<<c)))
            throw std::runtime_error("Nonfinite or duplicate cell");
        rows[r][c]=value;masks[r]|=1u<<c;
    };
    for(size_t p=sheet;p+4<=book.size();){
        unsigned id=u16(book,p),len=u16(book,p+2);size_t d=p+4;
        if(d+len>book.size()) throw std::runtime_error("Truncated sheet record");
        if(id==0x203){
            if(len!=14) throw std::runtime_error("Invalid NUMBER record");
            double v;std::memcpy(&v,book.data()+d+6,8);put(u16(book,d),u16(book,d+2),v);
        }else if(id==0x27e){
            if(len!=10) throw std::runtime_error("Invalid RK record");
            put(u16(book,d),u16(book,d+2),rk(u32(book,d+6)));
        }else if(id==0xbd){
            unsigned r=u16(book,d),first=u16(book,d+2),last=u16(book,d+len-2);
            if(last<first || len!=6+6*(last-first+1)) throw std::runtime_error("Invalid MULRK");
            for(unsigned c=first;c<=last;++c) put(r,c,rk(u32(book,d+6+6*(c-first))));
        }else if(id==6) throw std::runtime_error("Formula cells unsupported");
        else if(id==0x0a) break;
        p=d+len;
    }
    InputData result;
    result.dates.reserve(rows.size()-1);result.values.reserve(rows.size()-1);
    for(size_t r=1;r<rows.size();++r){const auto&v=rows[r];
        if(masks[r]!=511 || v[0]!=std::floor(v[0]) || v[0]<25569 || v[0]>100000)
            throw std::runtime_error("Missing row/cell or invalid date");
        if(!result.dates.empty() && v[0]<=result.dates.back())
            throw std::runtime_error("Dates must be strictly increasing");
        result.dates.push_back(int(v[0]));
        std::array<double,8> a;for(int j=0;j<8;++j)a[j]=v[j+1];
        result.values.push_back(a);
    }
    if(result.values.size()<2) throw std::runtime_error("Insufficient observations");
    return result;
}

