// Native boundary tests for the SAME generated frontend/header; no codec calls.
#define main MetadataDriverMain
#include "paired_metadata.cpp"
#undef main

template<class F> void Refuses(F action)
{
    bool refused = false;
    try { action(); } catch (const std::runtime_error&) { refused = true; }
    Require(refused, "metadata boundary did not refuse");
}

int main()
{
    try {
        namespace m = paired_metadata;
        unsigned checks = 0;
        Require(m::End(UINT64_MAX-1,1)==UINT64_MAX,"valid end boundary"); ++checks;
        Refuses([] { m::End(UINT64_MAX,1); }); ++checks;
        m::records = m::Store{};
        Refuses([] { m::Select(1,1,0,3,0,0); }); ++checks;
        Refuses([] { m::Select(0,1,0,3,0,0); }); ++checks;
        for (unsigned i=0;i<104;++i) m::Select(i,i%4,0,3,i%2,0);
        Require(m::records.selection_count==104,"selection capacity accepted"); ++checks;
        Refuses([] { m::Select(104,0,0,3,0,0); }); ++checks;
        Require(m::records.selection_count==104,"capacity failure is atomic"); ++checks;
        Buffer source(64), reference(64), scratch(64);
        std::vector<const void*> in(4097,source.data);
        std::vector<void*> out(1,scratch.data);
        m::records.snapshot_count=0;
        Refuses([&] { m::Capture(0,source,reference,scratch,NULL,scratch.data,64,in,out,64); }); ++checks;
        in.resize(1); out.resize(1025,scratch.data);
        Refuses([&] { m::Capture(0,source,reference,scratch,NULL,scratch.data,64,in,out,64); }); ++checks;
        out.resize(1);
        Refuses([&] { m::Capture(2,source,reference,scratch,NULL,scratch.data,64,in,out,64); }); ++checks;
        in.clear();
        Refuses([&] { m::Capture(0,source,reference,scratch,NULL,scratch.data,64,in,out,64); }); ++checks;
        Require(m::records.snapshot_count==0,"early capture failure is atomic"); ++checks;
        ElfW(Phdr) loads[17] = {};
        for (unsigned i=0;i<17;++i) {
            loads[i].p_type=PT_LOAD; loads[i].p_vaddr=0x1000; loads[i].p_memsz=0x100; loads[i].p_flags=PF_X;
        }
        dl_phdr_info info = {}; info.dlpi_phdr=loads; info.dlpi_phnum=16;
        m::Snapshot snapshot = {}; snapshot.main_address=0x1000;
        m::ImageQuery query = {&snapshot,false};
        Require(m::FindImage(&info,sizeof(info),&query)==1 && snapshot.segment_count==16,"segment capacity"); ++checks;
        Refuses([&] { m::FindImage(&info,sizeof(info),&query); }); ++checks;
        snapshot.segment_count=0; query.found=false; info.dlpi_phnum=17;
        Refuses([&] { m::FindImage(&info,sizeof(info),&query); }); ++checks;
        snapshot.segment_count=0; query.found=false; info.dlpi_phnum=1; info.dlpi_addr=UINT64_MAX-8;
        Refuses([&] { m::FindImage(&info,sizeof(info),&query); }); ++checks;
        m::Snapshot a = {}, b = {};
        Require(m::Same(a,b),"equal empty snapshot"); ++checks;
        a.input_count=b.input_count=1; b.inputs[0]=1;
        Require(!m::Same(a,b),"input pointer change"); ++checks; b.inputs[0]=0;
        a.output_count=b.output_count=1; b.outputs[0]=1;
        Require(!m::Same(a,b),"output pointer change"); ++checks; b.outputs[0]=0;
        a.segment_count=b.segment_count=1; b.segments[0].flags=1;
        Require(!m::Same(a,b),"segment change"); ++checks; b.segments[0].flags=0;
        b.scratch.bytes=1; Require(!m::Same(a,b),"allocation change"); ++checks; b.scratch.bytes=0;
        b.main_address=1; Require(!m::Same(a,b),"function change"); ++checks;
        std::printf("{\"schema\":\"paired-metadata-unit/v1\",\"cases\":%u,\"timed\":false}\n",checks);
        return 0;
    } catch (const std::exception& error) { std::fprintf(stderr,"%s\n",error.what()); return 1; }
}
