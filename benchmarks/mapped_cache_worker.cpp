// Same sequence-only workload compiled against v1.0.0 and the current library.
#include "Cache.h"
#include <htslib/faidx.h>
#include <openssl/evp.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point a) {return std::chrono::duration<double,std::milli>(Clock::now()-a).count();}
int main(int argc,char** argv) try {
    if(argc!=7)throw std::runtime_error("build|query|hold cache|indexed input output_or_ids count length");
    std::string mode=argv[1],arm=argv[2]; size_t count=std::stoull(argv[5]),length=std::stoull(argv[6]);
    TracEon::Cache cache;
    std::unique_ptr<faidx_t,decltype(&fai_destroy)> index(nullptr,fai_destroy);
    auto begin=Clock::now();
    if(mode=="build") {
        if(arm=="indexed") {if(fai_build(argv[3]))throw std::runtime_error("fai build");}
        else cache.loadFile(argv[3]);
        double load=elapsed(begin); begin=Clock::now();
        if(arm=="cache") {if(cache.size()!=count)throw std::runtime_error("count");cache.save(argv[4]);}
        std::cout<<"{\"load_ms\":"<<load<<",\"save_ms\":"<<elapsed(begin)<<"}\n";return 0;
    }
    if(arm=="indexed") {index.reset(fai_load3(argv[3],nullptr,nullptr,0));if(!index)throw std::runtime_error("index load");}
    else {cache.restore(argv[3]);if(cache.size()!=count)throw std::runtime_error("count");}
    double setup=elapsed(begin);
    std::unique_ptr<char,decltype(&free)> buffer(nullptr,free);
    auto get=[&](const std::string& id)->std::string_view {
        std::string_view s;
        if(arm=="indexed") {hts_pos_t n=0;buffer.reset(faidx_fetch_seq64(index.get(),id.c_str(),0,HTS_POS_MAX,&n));if(!buffer||n<0)throw std::runtime_error("miss");s={buffer.get(),size_t(n)};}
        else s=cache.getView(id);
        if(s.size()!=length)throw std::runtime_error("length/miss");return s;
    };
    std::vector<std::string> ids;
    if(mode=="query") {std::ifstream f(argv[4]);std::string s;while(std::getline(f,s))ids.push_back(s);if(!f.eof()||ids.size()<10000)throw std::runtime_error("ids");}
    double lookup=0;size_t sum=0;
    if(mode=="query" && arm=="cache") {
        for(const auto& id:ids)get(id); // Match the historical warmed lookup experiment.
        begin=Clock::now();for(int j=0;j<5;++j)for(const auto& id:ids)sum+=get(id).size();lookup=elapsed(begin);
    }
    std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> h(EVP_MD_CTX_new(),EVP_MD_CTX_free);
    if(!h||EVP_DigestInit_ex(h.get(),EVP_sha256(),nullptr)!=1)throw std::runtime_error("SHA init");
    auto field=[&](std::string_view s) {unsigned char n[8];for(int j=0;j<8;++j)n[j]=uint64_t(s.size())>>(8*j);if(EVP_DigestUpdate(h.get(),n,8)!=1||EVP_DigestUpdate(h.get(),s.data(),s.size())!=1)throw std::runtime_error("SHA update");};
    begin=Clock::now();size_t queries=mode=="hold"?count:10000;
    for(size_t i=0;i<queries;++i) {std::string id=mode=="hold"?"read"+std::to_string(i):ids[i];auto s=get(id);field(id);field(s);}
    double consume=elapsed(begin);unsigned char digest[32];unsigned n;
    if(EVP_DigestFinal_ex(h.get(),digest,&n)!=1||n!=32)throw std::runtime_error("SHA final");
    std::cout<<std::setprecision(17)<<"{\"setup_ms\":"<<setup<<",\"lookup_ms\":"<<lookup<<",\"lookup_queries\":"<<5*ids.size()<<",\"length_sum\":"<<sum<<",\"consume_ms\":"<<consume<<",\"records\":"<<queries<<",\"sha256\":\"";
    for(auto b:digest)std::cout<<std::hex<<std::setw(2)<<std::setfill('0')<<unsigned(b);
    std::cout<<"\"}"<<std::endl;
    if(mode=="hold") {std::string s;std::getline(std::cin,s);}
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
