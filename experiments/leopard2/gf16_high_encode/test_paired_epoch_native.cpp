// Same generated header/frontend, no public codec calls. Epoch boundary tests.
#define main MetadataDriverMain
#include "paired_metadata.cpp"
#undef main

template<class F> void Refuses(F action)
{
    bool refused = false;
    try { action(); } catch (const std::runtime_error&) { refused = true; }
    Require(refused, "epoch boundary did not refuse");
}

int main(int argc, char** argv)
{
    try {
        if (argc == 2) {
            const std::string mode(argv[1]);
            if (mode == "--mark-order") { LeoPairedWitnessMark(1); return 99; }
            if (mode == "--mark-duplicate") {
                LeoPairedWitnessMark(0); LeoPairedWitnessMark(0); return 99;
            }
            Require(mode == "--marks" || mode == "--mark-capacity", "unit mode");
            for (unsigned i=0; i<6; ++i) LeoPairedWitnessMark(i);
            if (mode == "--mark-capacity") { LeoPairedWitnessMark(6); return 99; }
            std::puts("{\"schema\":\"paired-epoch-unit/v1\",\"cases\":1,\"timed\":false}");
            return 0;
        }
        Require(argc == 1, "unit arguments");
        namespace m = paired_metadata;
        unsigned checks = 0;
        Require(m::End(UINT64_MAX-1,1)==UINT64_MAX,"valid end boundary"); ++checks;
        Refuses([] { m::End(UINT64_MAX,1); }); ++checks;
        m::records = m::Store{};
        Refuses([] { m::Select(1,1,0,3,0,0); }); ++checks;
        Refuses([] { m::Select(0,1,0,3,0,0); }); ++checks;
        for (unsigned i=0;i<312;++i) m::Select(i,i%4,0,3,i%2,0);
        Require(m::records.selection_count==312,"selection capacity accepted"); ++checks;
        Refuses([] { m::Select(312,0,0,3,0,0); }); ++checks;
        Require(m::records.selection_count==312,"capacity failure is atomic"); ++checks;
        Buffer source(64), reference(64), scratch(64);
        std::vector<const void*> in(4097,source.data);
        std::vector<void*> out(1,scratch.data);
        const auto capture = [&](unsigned endpoint) {
            m::Capture(endpoint,source,reference,scratch,NULL,scratch.data,64,in,out,64);
        };
        Refuses([&] { capture(0); }); ++checks;
        in.resize(1); out.resize(1025,scratch.data);
        Refuses([&] { capture(0); }); ++checks;
        out.resize(1);
        Refuses([&] { capture(2); }); ++checks;
        in.clear();
        Refuses([&] { capture(0); }); ++checks;
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
        in.resize(1,source.data);
        for (unsigned endpoint=0;endpoint<6;++endpoint) capture(endpoint);
        Require(m::records.snapshot_count==6,"six snapshots accepted"); ++checks;
        Refuses([&] { capture(6); }); ++checks;
        Require(m::records.snapshot_count==6,"snapshot capacity refusal is atomic"); ++checks;
        // A different but internally valid buffer in a later epoch must fail
        // against snapshot0, even when that epoch's before/after would agree.
        Buffer moved(64);
        for (unsigned endpoint : {2U,3U,4U,5U}) {
            m::records = m::Store{};
            in[0]=source.data;
            for (unsigned i=0;i<endpoint;++i) capture(i);
            in[0]=moved.data;
            Refuses([&] { m::Capture(endpoint,moved,reference,scratch,NULL,scratch.data,64,in,out,64); }); ++checks;
            // Captured != validated: the failed equality check follows capture.
            Require(m::records.snapshot_count==endpoint+1,"failed snapshot progress"); ++checks;
        }
#ifndef LEO_PAIRED_NATIVE
        namespace d = leopard2_internal;
        for (unsigned epoch=0;epoch<3;++epoch) {
            Require(d::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true),"arm late probe");
            // The old predicate really passes with the probe armed. No codec
            // calls are required to expose the zero-count final-slot hole.
            Require(d::AutoGF16GFNIEncodeCallCountForDiagnostics()==0 &&
                d::AutoGF16GFNIEncodeModeForDiagnostics()==1,"old gate false pass"); ++checks;
            Refuses([&] { m::RequireQuiescentProbe(0); }); ++checks;
            m::RequireQuiescentProbe(0); m::RequireQuiescentProbe(0);
            Require(!d::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(),"normal guard does not arm probe"); ++checks;
        }
        Require(d::SetAutoGF16GFNIEncodeEnabledForDiagnostics(false),"disabled probe");
        Refuses([&] { m::RequireQuiescentProbe(0); }); ++checks;
        Require(d::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(),"normalize disabled probe");
        Refuses([&] { m::RequireQuiescentProbe(0); }); ++checks;
        Require(d::SetAutoGF16GFNIEncodeEnabledForDiagnostics(true) &&
            d::FinishAutoGF16GFNIEncodeRouteProbeForDiagnostics(),"restore normal probe state");
        Refuses([&] { m::RequireQuiescentProbe(1); }); ++checks;
#endif
        m::records = m::Store{}; // Unit manipulations must not impersonate run progress.
        std::printf("{\"schema\":\"paired-epoch-unit/v1\",\"cases\":%u,\"timed\":false}\n",checks);
        return 0;
    } catch (const std::exception& error) { std::fprintf(stderr,"%s\n",error.what()); return 1; }
}
