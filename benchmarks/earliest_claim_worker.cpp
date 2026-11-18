// Historical Cache linked from fa95f1d, not the current implementation.
#include "Cache.h"
#include <htslib/faidx.h>
#include <openssl/evp.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <vector>
#include <cstdlib>
// Reuse the historical scanner verbatim, retaining its full-file traversal.
#define main historical_disk_main
#include "benchmark_disk.cpp"
#undef main
using Clock=std::chrono::steady_clock;
static double ms(Clock::time_point a,Clock::time_point b) {return std::chrono::duration<double,std::milli>(b-a).count();}
static void update(EVP_MD_CTX* h,const void* p,size_t n) {if(EVP_DigestUpdate(h,p,n)!=1)throw std::runtime_error("digest update");}
static void field(EVP_MD_CTX* h,const std::string& s) {unsigned char n[8];for(int i=0;i<8;++i)n[i]=uint64_t(s.size())>>(i*8);update(h,n,8);update(h,s.data(),s.size());}
int main(int argc,char** argv) try {
    if(argc!=5)throw std::runtime_error("worker scan|cache|indexed|roundtrip reads ids cache_or_fai");
    std::string arm=argv[1];
    if(arm!="scan"&&arm!="cache"&&arm!="indexed"&&arm!="roundtrip")throw std::runtime_error("arm");
    std::vector<std::string> ids;std::ifstream in(argv[3]);std::string id;
    while(std::getline(in,id)) {if(id.empty())throw std::runtime_error("empty ID");ids.push_back(id);}
    if(!in.eof()||ids.empty())throw std::runtime_error("ID file");
    TracEon::Cache cache;
    std::unique_ptr<faidx_t,decltype(&fai_destroy)> index(nullptr,fai_destroy);
    auto begin=Clock::now();
    if(arm=="cache"||arm=="roundtrip")cache.loadFasta(argv[2]);
    if(arm=="indexed") {index.reset(fai_load3(argv[2],argv[4],nullptr,0));if(!index)throw std::runtime_error("index");}
    auto ready=Clock::now();double save_ms=0,restore_ms=0;
    if(arm=="roundtrip") {
        auto a=Clock::now();cache.save(argv[4]);auto b=Clock::now();
        cache.restore(argv[4]);auto c=Clock::now();save_ms=ms(a,b);restore_ms=ms(b,c);
    }
    std::unique_ptr<EVP_MD_CTX,decltype(&EVP_MD_CTX_free)> h(EVP_MD_CTX_new(),EVP_MD_CTX_free);
    if(!h||EVP_DigestInit_ex(h.get(),EVP_sha256(),nullptr)!=1)throw std::runtime_error("digest init");
    size_t bytes=0;auto start=Clock::now();
    for(const auto& key:ids) {
        std::string seq;
        if(arm=="scan")seq=find_sequence_on_disk(argv[2],key);
        else if(arm=="indexed") {
            hts_pos_t n=0;std::unique_ptr<char,decltype(&free)> p(faidx_fetch_seq64(index.get(),key.c_str(),0,HTS_POS_MAX,&n),free);
            if(!p||n<0)throw std::runtime_error("indexed miss");seq.assign(p.get(),n);
        } else seq=cache.get(key);
        if(seq.empty())throw std::runtime_error("empty result or miss");
        field(h.get(),key);field(h.get(),seq);bytes+=seq.size();
    }
    auto end=Clock::now();unsigned char digest[32];unsigned n;
    if(EVP_DigestFinal_ex(h.get(),digest,&n)!=1||n!=32)throw std::runtime_error("digest final");
    std::cout<<std::setprecision(17)<<"{\"setup_ms\":"<<ms(begin,ready)<<",\"query_ms\":"<<ms(start,end)
        <<",\"total_ms\":"<<ms(begin,ready)+ms(start,end)<<",\"save_ms\":"<<save_ms<<",\"restore_ms\":"<<restore_ms
        <<",\"records\":"<<ids.size()<<",\"cache_count\":"<<cache.size()<<",\"bytes\":"<<bytes<<",\"sha256\":\"";
    for(auto b:digest)std::cout<<std::hex<<std::setfill('0')<<std::setw(2)<<unsigned(b);
    std::cout<<"\"}\n";
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
