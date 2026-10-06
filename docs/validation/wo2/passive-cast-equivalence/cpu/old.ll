; ModuleID = '/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/passive-cast-proof/probe.cpp'
source_filename = "/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/passive-cast-proof/probe.cpp"
target datalayout = "e-m:e-p270:32:32-p271:32:32-p272:64:64-i64:64-i128:128-f80:128-n8:16:32:64-S128"
target triple = "x86_64-unknown-linux-gnu"

%"class.std::map" = type { %"class.std::_Rb_tree" }
%"class.std::_Rb_tree" = type { %"struct.std::_Rb_tree<std::__cxx11::basic_string<char>, std::pair<const std::__cxx11::basic_string<char>, Kokkos::Tools::Experimental::TeamSizeTuner>, std::_Select1st<std::pair<const std::__cxx11::basic_string<char>, Kokkos::Tools::Experimental::TeamSizeTuner>>, std::less<std::__cxx11::basic_string<char>>>::_Rb_tree_impl" }
%"struct.std::_Rb_tree<std::__cxx11::basic_string<char>, std::pair<const std::__cxx11::basic_string<char>, Kokkos::Tools::Experimental::TeamSizeTuner>, std::_Select1st<std::pair<const std::__cxx11::basic_string<char>, Kokkos::Tools::Experimental::TeamSizeTuner>>, std::less<std::__cxx11::basic_string<char>>>::_Rb_tree_impl" = type { [8 x i8], %"struct.std::_Rb_tree_header" }
%"struct.std::_Rb_tree_header" = type { %"struct.std::_Rb_tree_node_base", i64 }
%"struct.std::_Rb_tree_node_base" = type { i32, ptr, ptr, ptr }

$_ZNSt3mapINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEN6Kokkos5Tools12Experimental13TeamSizeTunerESt4lessIS5_ESaISt4pairIKS5_S9_EEED2Ev = comdat any

$__clang_call_terminate = comdat any

$_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE8_M_eraseEPSt13_Rb_tree_nodeISC_E = comdat any

@_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E = internal global %"class.std::map" zeroinitializer, align 8
@__dso_handle = external hidden global i8
@llvm.global_ctors = appending global [1 x { i32, ptr, ptr }] [{ i32, ptr, ptr } { i32 65535, ptr @_GLOBAL__sub_I_probe.cpp, ptr null }]
@llvm.compiler.used = appending global [2 x ptr] [ptr @wo2_passive_encode, ptr @wo2_passive_specific], section "llvm.metadata"

; Function Attrs: mustprogress nounwind uwtable
define linkonce_odr dso_local void @_ZNSt3mapINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEN6Kokkos5Tools12Experimental13TeamSizeTunerESt4lessIS5_ESaISt4pairIKS5_S9_EEED2Ev(ptr noundef nonnull align 8 dereferenceable(48) %this) unnamed_addr #0 comdat align 2 personality ptr @__gxx_personality_v0 !dbg !6 {
entry:
  %_M_parent.i.i.i = getelementptr inbounds nuw i8, ptr %this, i64 16, !dbg !10 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:733:64
  %0 = load ptr, ptr %_M_parent.i.i.i, align 8, !dbg !10, !tbaa !18 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:733:64
  invoke void @_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE8_M_eraseEPSt13_Rb_tree_nodeISC_E(ptr noundef nonnull align 8 dereferenceable(48) %this, ptr noundef %0)
          to label %_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EED2Ev.exit unwind label %terminate.lpad.i, !dbg !27 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:982:9

terminate.lpad.i:                                 ; preds = %entry
  %1 = landingpad { ptr, i32 }
          catch ptr null, !dbg !27 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:982:9
  %2 = extractvalue { ptr, i32 } %1, 0, !dbg !27 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:982:9
  tail call void @__clang_call_terminate(ptr %2) #11, !dbg !27 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:982:9
  unreachable, !dbg !27 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:982:9

_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EED2Ev.exit: ; preds = %entry
  ret void, !dbg !28 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_map.h:314:22
}

; Function Attrs: nofree nounwind
declare dso_local i32 @__cxa_atexit(ptr, ptr, ptr) local_unnamed_addr #1

; Function Attrs: mustprogress nofree noinline nounwind willreturn memory(write) uwtable
define dso_local void @wo2_passive_specific(double noundef %rho, double noundef %ppar, double noundef %pperp, double noundef %bmag, ptr nocapture noundef writeonly initializes((0, 16)) %out) #2 !dbg !29 {
entry:
  %call.i = tail call double @log(double noundef %rho) #12, !dbg !31, !tbaa !35 ; include/eos/cgl_passive.hpp:18:21
  %call1.i = tail call double @log(double noundef %bmag) #12, !dbg !37, !tbaa !35 ; include/eos/cgl_passive.hpp:18:38
  %call2.i = tail call double @log(double noundef %ppar) #12, !dbg !38, !tbaa !35 ; include/eos/cgl_passive.hpp:19:11
  %0 = tail call double @llvm.fmuladd.f64(double %call1.i, double 2.000000e+00, double %call2.i), !dbg !39 ; include/eos/cgl_passive.hpp:19:21
  %1 = tail call double @llvm.fmuladd.f64(double %call.i, double -3.000000e+00, double %0), !dbg !40 ; include/eos/cgl_passive.hpp:19:32
  %call3.i = tail call double @log(double noundef %pperp) #12, !dbg !41, !tbaa !35 ; include/eos/cgl_passive.hpp:20:11
  %call4.i = tail call double @log(double noundef %ppar) #12, !dbg !42, !tbaa !35 ; include/eos/cgl_passive.hpp:20:24
  %sub.i = fsub double %call3.i, %call4.i, !dbg !43 ; include/eos/cgl_passive.hpp:20:22
  %2 = tail call double @llvm.fmuladd.f64(double %call.i, double 2.000000e+00, double %sub.i), !dbg !44 ; include/eos/cgl_passive.hpp:20:34
  %3 = tail call double @llvm.fmuladd.f64(double %call1.i, double -3.000000e+00, double %2), !dbg !45 ; include/eos/cgl_passive.hpp:20:45
  store double %1, ptr %out, align 8, !dbg !46, !tbaa !47 ; probe.cpp:13:10
  %arrayidx1 = getelementptr inbounds nuw i8, ptr %out, i64 8, !dbg !49 ; probe.cpp:14:3
  store double %3, ptr %arrayidx1, align 8, !dbg !50, !tbaa !47 ; probe.cpp:14:10
  ret void, !dbg !51 ; probe.cpp:15:1
}

; Function Attrs: mustprogress nofree noinline nounwind willreturn memory(write) uwtable
define dso_local void @wo2_passive_encode(double noundef %rho, double noundef %ppar, double noundef %pperp, double noundef %bmag, ptr nocapture noundef writeonly initializes((0, 16)) %out) #2 !dbg !52 {
entry:
  %call.i.i = tail call double @log(double noundef %rho) #12, !dbg !53, !tbaa !35 ; include/eos/cgl_passive.hpp:18:21
  %call1.i.i = tail call double @log(double noundef %bmag) #12, !dbg !57, !tbaa !35 ; include/eos/cgl_passive.hpp:18:38
  %call2.i.i = tail call double @log(double noundef %ppar) #12, !dbg !58, !tbaa !35 ; include/eos/cgl_passive.hpp:19:11
  %0 = tail call double @llvm.fmuladd.f64(double %call1.i.i, double 2.000000e+00, double %call2.i.i), !dbg !59 ; include/eos/cgl_passive.hpp:19:21
  %1 = tail call double @llvm.fmuladd.f64(double %call.i.i, double -3.000000e+00, double %0), !dbg !60 ; include/eos/cgl_passive.hpp:19:32
  %call3.i.i = tail call double @log(double noundef %pperp) #12, !dbg !61, !tbaa !35 ; include/eos/cgl_passive.hpp:20:11
  %call4.i.i = tail call double @log(double noundef %ppar) #12, !dbg !62, !tbaa !35 ; include/eos/cgl_passive.hpp:20:24
  %sub.i.i = fsub double %call3.i.i, %call4.i.i, !dbg !63 ; include/eos/cgl_passive.hpp:20:22
  %2 = tail call double @llvm.fmuladd.f64(double %call.i.i, double 2.000000e+00, double %sub.i.i), !dbg !64 ; include/eos/cgl_passive.hpp:20:34
  %3 = tail call double @llvm.fmuladd.f64(double %call1.i.i, double -3.000000e+00, double %2), !dbg !65 ; include/eos/cgl_passive.hpp:20:45
  %mul.i = fmul double %rho, %1, !dbg !66 ; include/eos/cgl_passive.hpp:27:14
  %mul3.i = fmul double %rho, %3, !dbg !67 ; include/eos/cgl_passive.hpp:27:23
  store double %mul.i, ptr %out, align 8, !dbg !68, !tbaa !47 ; probe.cpp:19:10
  %arrayidx1 = getelementptr inbounds nuw i8, ptr %out, i64 8, !dbg !69 ; probe.cpp:20:3
  store double %mul3.i, ptr %arrayidx1, align 8, !dbg !70, !tbaa !47 ; probe.cpp:20:10
  ret void, !dbg !71 ; probe.cpp:21:1
}

declare dso_local i32 @__gxx_personality_v0(...)

; Function Attrs: noinline noreturn nounwind uwtable
define linkonce_odr hidden void @__clang_call_terminate(ptr noundef %0) local_unnamed_addr #3 comdat {
  %2 = tail call ptr @__cxa_begin_catch(ptr %0) #12
  tail call void @_ZSt9terminatev() #11
  unreachable
}

declare dso_local ptr @__cxa_begin_catch(ptr) local_unnamed_addr

; Function Attrs: cold nofree noreturn
declare dso_local void @_ZSt9terminatev() local_unnamed_addr #4

; Function Attrs: mustprogress uwtable
define linkonce_odr dso_local void @_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE8_M_eraseEPSt13_Rb_tree_nodeISC_E(ptr noundef nonnull align 8 dereferenceable(48) %this, ptr noundef %__x) local_unnamed_addr #5 comdat align 2 personality ptr @__gxx_personality_v0 !dbg !72 {
entry:
  %cmp.not6 = icmp eq ptr %__x, null, !dbg !73, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:1930:18 at depth 1
  br i1 %cmp.not6, label %while.end, label %while.body, !dbg !75, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:1930:7 at depth 1

while.body:                                       ; preds = %entry, %_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE12_M_drop_nodeEPSt13_Rb_tree_nodeISC_E.exit
  %__x.addr.07 = phi ptr [ %1, %_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE12_M_drop_nodeEPSt13_Rb_tree_nodeISC_E.exit ], [ %__x, %entry ]
  %_M_right.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 24, !dbg !76 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:786:45
  %0 = load ptr, ptr %_M_right.i, align 8, !dbg !76, !tbaa !79 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:786:45
  tail call void @_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE8_M_eraseEPSt13_Rb_tree_nodeISC_E(ptr noundef nonnull align 8 dereferenceable(48) %this, ptr noundef %0), !dbg !80, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:1932:4 at depth 1
  %_M_left.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 16, !dbg !81 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:778:45
  %1 = load ptr, ptr %_M_left.i, align 8, !dbg !81, !tbaa !84 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:778:45
  %_M_storage.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 32, !dbg !85 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:231:16
  %m_variable_names.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 136, !dbg !92 ; source-unfused/kokkos/core/src/Kokkos_Tuners.hpp:279:7
  %2 = load ptr, ptr %m_variable_names.i.i.i.i.i, align 8, !dbg !107, !tbaa !111 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:735:30
  %_M_finish.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 144, !dbg !114 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:735:54
  %3 = load ptr, ptr %_M_finish.i.i.i.i.i.i, align 8, !dbg !114, !tbaa !115 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:735:54
  %cmp.not3.i.i.i.i.i.i.i.i = icmp eq ptr %2, %3, !dbg !116, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:19
  br i1 %cmp.not3.i.i.i.i.i.i.i.i, label %invoke.cont.i.i.i.i.i.i, label %for.body.i.i.i.i.i.i.i.i, !dbg !124, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:4

for.body.i.i.i.i.i.i.i.i:                         ; preds = %while.body, %_ZSt8_DestroyINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEvPT_.exit.i.i.i.i.i.i.i.i
  %__first.addr.04.i.i.i.i.i.i.i.i = phi ptr [ %incdec.ptr.i.i.i.i.i.i.i.i, %_ZSt8_DestroyINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEvPT_.exit.i.i.i.i.i.i.i.i ], [ %2, %while.body ]
  %4 = load ptr, ptr %__first.addr.04.i.i.i.i.i.i.i.i, align 8, !dbg !125, !tbaa !137 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:228:28
  %5 = getelementptr inbounds nuw i8, ptr %__first.addr.04.i.i.i.i.i.i.i.i, i64 16, !dbg !141 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:246:57
  %cmp.i.i.i.i.i.i.i.i.i.i.i.i = icmp eq ptr %4, %5, !dbg !144 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:269:16
  br i1 %cmp.i.i.i.i.i.i.i.i.i.i.i.i, label %_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE11_M_is_localEv.exit.thread.i.i.i.i.i.i.i.i.i.i.i, label %if.then.i.i.i.i.i.i.i.i.i.i.i, !dbg !144 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:269:16

_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE11_M_is_localEv.exit.thread.i.i.i.i.i.i.i.i.i.i.i: ; preds = %for.body.i.i.i.i.i.i.i.i
  %_M_string_length.i.i.i.i.i.i.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__first.addr.04.i.i.i.i.i.i.i.i, i64 8, !dbg !145 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:271:10
  %6 = load i64, ptr %_M_string_length.i.i.i.i.i.i.i.i.i.i.i.i, align 8, !dbg !145, !tbaa !146 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:271:10
  %cmp3.i.i.i.i.i.i.i.i.i.i.i.i = icmp ult i64 %6, 16, !dbg !147 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:271:27
  tail call void @llvm.assume(i1 %cmp3.i.i.i.i.i.i.i.i.i.i.i.i), !dbg !147 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:271:27
  br label %_ZSt8_DestroyINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEvPT_.exit.i.i.i.i.i.i.i.i, !dbg !148 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:287:6

if.then.i.i.i.i.i.i.i.i.i.i.i:                    ; preds = %for.body.i.i.i.i.i.i.i.i
  %7 = load i64, ptr %5, align 8, !dbg !149, !tbaa !150 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:288:15
  %8 = add i64 %7, 1, !dbg !151 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:294:73
  tail call void @_ZdlPvm(ptr noundef %4, i64 noundef %8) #13, !dbg !154 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/new_allocator.h:172:2
  br label %_ZSt8_DestroyINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEvPT_.exit.i.i.i.i.i.i.i.i, !dbg !159 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:288:4

_ZSt8_DestroyINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEvPT_.exit.i.i.i.i.i.i.i.i: ; preds = %if.then.i.i.i.i.i.i.i.i.i.i.i, %_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE11_M_is_localEv.exit.thread.i.i.i.i.i.i.i.i.i.i.i
  %incdec.ptr.i.i.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__first.addr.04.i.i.i.i.i.i.i.i, i64 32, !dbg !160 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:30
  %cmp.not.i.i.i.i.i.i.i.i = icmp eq ptr %incdec.ptr.i.i.i.i.i.i.i.i, %3, !dbg !116, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:19
  br i1 %cmp.not.i.i.i.i.i.i.i.i, label %invoke.contthread-pre-split.i.i.i.i.i.i, label %for.body.i.i.i.i.i.i.i.i, !dbg !124, !llvm.loop !161, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:4

invoke.contthread-pre-split.i.i.i.i.i.i:          ; preds = %_ZSt8_DestroyINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEEvPT_.exit.i.i.i.i.i.i.i.i
  %.pr.i.i.i.i.i.i = load ptr, ptr %m_variable_names.i.i.i.i.i, align 8, !dbg !164, !tbaa !111 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:368:24
  br label %invoke.cont.i.i.i.i.i.i, !dbg !164 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:368:24

invoke.cont.i.i.i.i.i.i:                          ; preds = %invoke.contthread-pre-split.i.i.i.i.i.i, %while.body
  %9 = phi ptr [ %.pr.i.i.i.i.i.i, %invoke.contthread-pre-split.i.i.i.i.i.i ], [ %2, %while.body ], !dbg !164 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:368:24
  %tobool.not.i.i.i.i.i.i.i.i = icmp eq ptr %9, null, !dbg !167 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:388:6
  br i1 %tobool.not.i.i.i.i.i.i.i.i, label %_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev.exit.i.i.i.i.i, label %if.then.i.i.i.i.i.i.i.i, !dbg !167 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:388:6

if.then.i.i.i.i.i.i.i.i:                          ; preds = %invoke.cont.i.i.i.i.i.i
  %_M_end_of_storage.i.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 152, !dbg !170 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:17
  %10 = load ptr, ptr %_M_end_of_storage.i.i.i.i.i.i.i, align 8, !dbg !170, !tbaa !171 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:17
  %sub.ptr.lhs.cast.i.i.i.i.i.i.i = ptrtoint ptr %10 to i64, !dbg !172 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  %sub.ptr.rhs.cast.i.i.i.i.i.i.i = ptrtoint ptr %9 to i64, !dbg !172 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  %sub.ptr.sub.i.i.i.i.i.i.i = sub i64 %sub.ptr.lhs.cast.i.i.i.i.i.i.i, %sub.ptr.rhs.cast.i.i.i.i.i.i.i, !dbg !172 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  tail call void @_ZdlPvm(ptr noundef nonnull %9, i64 noundef %sub.ptr.sub.i.i.i.i.i.i.i) #13, !dbg !173 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/new_allocator.h:172:2
  br label %_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev.exit.i.i.i.i.i, !dbg !178 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:389:4

_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev.exit.i.i.i.i.i: ; preds = %if.then.i.i.i.i.i.i.i.i, %invoke.cont.i.i.i.i.i.i
  %m_space.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 72, !dbg !92 ; source-unfused/kokkos/core/src/Kokkos_Tuners.hpp:279:7
  %sub_values.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 96, !dbg !179 ; source-unfused/kokkos/core/src/Kokkos_Tuners.hpp:69:8
  %11 = load ptr, ptr %sub_values.i.i.i.i.i.i, align 8, !dbg !182, !tbaa !185 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:735:30
  %_M_finish.i.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 104, !dbg !188 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:735:54
  %12 = load ptr, ptr %_M_finish.i.i.i.i.i.i.i, align 8, !dbg !188, !tbaa !189 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:735:54
  %cmp.not3.i.i.i.i.i.i.i.i.i = icmp eq ptr %11, %12, !dbg !190, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:19
  br i1 %cmp.not3.i.i.i.i.i.i.i.i.i, label %invoke.cont.i.i.i.i.i.i.i, label %for.body.i.i.i.i.i.i.i.i.i, !dbg !197, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:4

for.body.i.i.i.i.i.i.i.i.i:                       ; preds = %_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev.exit.i.i.i.i.i, %_ZSt8_DestroyIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEEEvPT_.exit.i.i.i.i.i.i.i.i.i
  %__first.addr.04.i.i.i.i.i.i.i.i.i = phi ptr [ %incdec.ptr.i.i.i.i.i.i.i.i.i, %_ZSt8_DestroyIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEEEvPT_.exit.i.i.i.i.i.i.i.i.i ], [ %11, %_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev.exit.i.i.i.i.i ]
  %13 = load ptr, ptr %__first.addr.04.i.i.i.i.i.i.i.i.i, align 8, !dbg !198, !tbaa !207 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:368:24
  %tobool.not.i.i.i.i.i.i.i.i.i.i.i.i.i.i = icmp eq ptr %13, null, !dbg !210 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:388:6
  br i1 %tobool.not.i.i.i.i.i.i.i.i.i.i.i.i.i.i, label %_ZSt8_DestroyIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEEEvPT_.exit.i.i.i.i.i.i.i.i.i, label %if.then.i.i.i.i.i.i.i.i.i.i.i.i.i.i, !dbg !210 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:388:6

if.then.i.i.i.i.i.i.i.i.i.i.i.i.i.i:              ; preds = %for.body.i.i.i.i.i.i.i.i.i
  %_M_end_of_storage.i.i.i.i.i.i.i.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__first.addr.04.i.i.i.i.i.i.i.i.i, i64 16, !dbg !213 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:17
  %14 = load ptr, ptr %_M_end_of_storage.i.i.i.i.i.i.i.i.i.i.i.i.i, align 8, !dbg !213, !tbaa !214 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:17
  %sub.ptr.lhs.cast.i.i.i.i.i.i.i.i.i.i.i.i.i = ptrtoint ptr %14 to i64, !dbg !215 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  %sub.ptr.rhs.cast.i.i.i.i.i.i.i.i.i.i.i.i.i = ptrtoint ptr %13 to i64, !dbg !215 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  %sub.ptr.sub.i.i.i.i.i.i.i.i.i.i.i.i.i = sub i64 %sub.ptr.lhs.cast.i.i.i.i.i.i.i.i.i.i.i.i.i, %sub.ptr.rhs.cast.i.i.i.i.i.i.i.i.i.i.i.i.i, !dbg !215 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  tail call void @_ZdlPvm(ptr noundef nonnull %13, i64 noundef %sub.ptr.sub.i.i.i.i.i.i.i.i.i.i.i.i.i) #13, !dbg !216 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/new_allocator.h:172:2
  br label %_ZSt8_DestroyIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEEEvPT_.exit.i.i.i.i.i.i.i.i.i, !dbg !221 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:389:4

_ZSt8_DestroyIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEEEvPT_.exit.i.i.i.i.i.i.i.i.i: ; preds = %if.then.i.i.i.i.i.i.i.i.i.i.i.i.i.i, %for.body.i.i.i.i.i.i.i.i.i
  %incdec.ptr.i.i.i.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__first.addr.04.i.i.i.i.i.i.i.i.i, i64 24, !dbg !222 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:30
  %cmp.not.i.i.i.i.i.i.i.i.i = icmp eq ptr %incdec.ptr.i.i.i.i.i.i.i.i.i, %12, !dbg !190, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:19
  br i1 %cmp.not.i.i.i.i.i.i.i.i.i, label %invoke.contthread-pre-split.i.i.i.i.i.i.i, label %for.body.i.i.i.i.i.i.i.i.i, !dbg !197, !llvm.loop !223, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h:162:4

invoke.contthread-pre-split.i.i.i.i.i.i.i:        ; preds = %_ZSt8_DestroyIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEEEvPT_.exit.i.i.i.i.i.i.i.i.i
  %.pr.i.i.i.i.i.i.i = load ptr, ptr %sub_values.i.i.i.i.i.i, align 8, !dbg !225, !tbaa !185 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:368:24
  br label %invoke.cont.i.i.i.i.i.i.i, !dbg !225 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:368:24

invoke.cont.i.i.i.i.i.i.i:                        ; preds = %invoke.contthread-pre-split.i.i.i.i.i.i.i, %_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev.exit.i.i.i.i.i
  %15 = phi ptr [ %.pr.i.i.i.i.i.i.i, %invoke.contthread-pre-split.i.i.i.i.i.i.i ], [ %11, %_ZNSt6vectorINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EED2Ev.exit.i.i.i.i.i ], !dbg !225 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:368:24
  %tobool.not.i.i.i.i.i.i.i.i.i = icmp eq ptr %15, null, !dbg !228 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:388:6
  br i1 %tobool.not.i.i.i.i.i.i.i.i.i, label %_ZNSt6vectorIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEESaIS5_EED2Ev.exit.i.i.i.i.i.i, label %if.then.i.i.i.i.i.i.i.i.i, !dbg !228 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:388:6

if.then.i.i.i.i.i.i.i.i.i:                        ; preds = %invoke.cont.i.i.i.i.i.i.i
  %_M_end_of_storage.i.i.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 112, !dbg !231 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:17
  %16 = load ptr, ptr %_M_end_of_storage.i.i.i.i.i.i.i.i, align 8, !dbg !231, !tbaa !232 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:17
  %sub.ptr.lhs.cast.i.i.i.i.i.i.i.i = ptrtoint ptr %16 to i64, !dbg !233 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  %sub.ptr.rhs.cast.i.i.i.i.i.i.i.i = ptrtoint ptr %15 to i64, !dbg !233 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  %sub.ptr.sub.i.i.i.i.i.i.i.i = sub i64 %sub.ptr.lhs.cast.i.i.i.i.i.i.i.i, %sub.ptr.rhs.cast.i.i.i.i.i.i.i.i, !dbg !233 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  tail call void @_ZdlPvm(ptr noundef nonnull %15, i64 noundef %sub.ptr.sub.i.i.i.i.i.i.i.i) #13, !dbg !234 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/new_allocator.h:172:2
  br label %_ZNSt6vectorIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEESaIS5_EED2Ev.exit.i.i.i.i.i.i, !dbg !239 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:389:4

_ZNSt6vectorIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEESaIS5_EED2Ev.exit.i.i.i.i.i.i: ; preds = %if.then.i.i.i.i.i.i.i.i.i, %invoke.cont.i.i.i.i.i.i.i
  %17 = load ptr, ptr %m_space.i.i.i.i.i, align 8, !dbg !240, !tbaa !207 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:368:24
  %tobool.not.i.i.i3.i.i.i.i.i.i = icmp eq ptr %17, null, !dbg !243 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:388:6
  br i1 %tobool.not.i.i.i3.i.i.i.i.i.i, label %_ZN6Kokkos5Tools12Experimental13TeamSizeTunerD2Ev.exit.i.i.i, label %if.then.i.i.i4.i.i.i.i.i.i, !dbg !243 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:388:6

if.then.i.i.i4.i.i.i.i.i.i:                       ; preds = %_ZNSt6vectorIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEESaIS5_EED2Ev.exit.i.i.i.i.i.i
  %_M_end_of_storage.i.i5.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 88, !dbg !245 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:17
  %18 = load ptr, ptr %_M_end_of_storage.i.i5.i.i.i.i.i.i, align 8, !dbg !245, !tbaa !214 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:17
  %sub.ptr.lhs.cast.i.i6.i.i.i.i.i.i = ptrtoint ptr %18 to i64, !dbg !246 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  %sub.ptr.rhs.cast.i.i7.i.i.i.i.i.i = ptrtoint ptr %17 to i64, !dbg !246 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  %sub.ptr.sub.i.i8.i.i.i.i.i.i = sub i64 %sub.ptr.lhs.cast.i.i6.i.i.i.i.i.i, %sub.ptr.rhs.cast.i.i7.i.i.i.i.i.i, !dbg !246 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:369:35
  tail call void @_ZdlPvm(ptr noundef nonnull %17, i64 noundef %sub.ptr.sub.i.i8.i.i.i.i.i.i) #13, !dbg !247 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/new_allocator.h:172:2
  br label %_ZN6Kokkos5Tools12Experimental13TeamSizeTunerD2Ev.exit.i.i.i, !dbg !250 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h:389:4

_ZN6Kokkos5Tools12Experimental13TeamSizeTunerD2Ev.exit.i.i.i: ; preds = %if.then.i.i.i4.i.i.i.i.i.i, %_ZNSt6vectorIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEESaIS5_EED2Ev.exit.i.i.i.i.i.i
  %19 = load ptr, ptr %_M_storage.i.i.i, align 8, !dbg !251, !tbaa !137 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:228:28
  %20 = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 48, !dbg !256 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:246:57
  %cmp.i.i.i.i.i.i = icmp eq ptr %19, %20, !dbg !258 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:269:16
  br i1 %cmp.i.i.i.i.i.i, label %_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE11_M_is_localEv.exit.thread.i.i.i.i.i, label %if.then.i.i.i.i.i, !dbg !258 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:269:16

_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE11_M_is_localEv.exit.thread.i.i.i.i.i: ; preds = %_ZN6Kokkos5Tools12Experimental13TeamSizeTunerD2Ev.exit.i.i.i
  %_M_string_length.i.i.i.i.i.i = getelementptr inbounds nuw i8, ptr %__x.addr.07, i64 40, !dbg !259 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:271:10
  %21 = load i64, ptr %_M_string_length.i.i.i.i.i.i, align 8, !dbg !259, !tbaa !146 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:271:10
  %cmp3.i.i.i.i.i.i = icmp ult i64 %21, 16, !dbg !260 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:271:27
  tail call void @llvm.assume(i1 %cmp3.i.i.i.i.i.i), !dbg !260 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:271:27
  br label %_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE12_M_drop_nodeEPSt13_Rb_tree_nodeISC_E.exit, !dbg !261 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:287:6

if.then.i.i.i.i.i:                                ; preds = %_ZN6Kokkos5Tools12Experimental13TeamSizeTunerD2Ev.exit.i.i.i
  %22 = load i64, ptr %20, align 8, !dbg !262, !tbaa !150 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:288:15
  %23 = add i64 %22, 1, !dbg !263 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:294:73
  tail call void @_ZdlPvm(ptr noundef %19, i64 noundef %23) #13, !dbg !265 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/new_allocator.h:172:2
  br label %_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE12_M_drop_nodeEPSt13_Rb_tree_nodeISC_E.exit, !dbg !268 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h:288:4

_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE12_M_drop_nodeEPSt13_Rb_tree_nodeISC_E.exit: ; preds = %_ZNKSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE11_M_is_localEv.exit.thread.i.i.i.i.i, %if.then.i.i.i.i.i
  tail call void @_ZdlPvm(ptr noundef nonnull %__x.addr.07, i64 noundef 168) #13, !dbg !269 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/new_allocator.h:172:2
  %cmp.not = icmp eq ptr %1, null, !dbg !73, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:1930:18 at depth 1
  br i1 %cmp.not, label %while.end, label %while.body, !dbg !75, !llvm.loop !276, !cray.depth !74 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:1930:7 at depth 1

while.end:                                        ; preds = %_ZNSt8_Rb_treeINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESt4pairIKS5_N6Kokkos5Tools12Experimental13TeamSizeTunerEESt10_Select1stISC_ESt4lessIS5_ESaISC_EE12_M_drop_nodeEPSt13_Rb_tree_nodeISC_E.exit, %entry
  ret void, !dbg !278 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:1937:5
}

; Function Attrs: nobuiltin nounwind
declare dso_local void @_ZdlPvm(ptr noundef, i64 noundef) local_unnamed_addr #6

; Function Attrs: mustprogress nofree nounwind willreturn memory(write)
declare dso_local double @log(double noundef) local_unnamed_addr #7

; Function Attrs: mustprogress nocallback nofree nosync nounwind speculatable willreturn memory(none)
declare double @llvm.fmuladd.f64(double, double, double) #8

; Function Attrs: nofree nounwind uwtable
define internal void @_GLOBAL__sub_I_probe.cpp() #9 section ".text.startup" personality ptr @__gxx_personality_v0 !dbg !279 {
entry:
  store i32 0, ptr getelementptr inbounds nuw (i8, ptr @_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E, i64 8), align 8, !dbg !280, !tbaa !293 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:171:26
  store ptr null, ptr getelementptr inbounds nuw (i8, ptr @_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E, i64 16), align 8, !dbg !294, !tbaa !18 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:204:27
  store ptr getelementptr inbounds nuw (i8, ptr @_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E, i64 8), ptr getelementptr inbounds nuw (i8, ptr @_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E, i64 24), align 8, !dbg !297, !tbaa !298 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:205:25
  store ptr getelementptr inbounds nuw (i8, ptr @_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E, i64 8), ptr getelementptr inbounds nuw (i8, ptr @_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E, i64 32), align 8, !dbg !299, !tbaa !300 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:206:26
  store i64 0, ptr getelementptr inbounds nuw (i8, ptr @_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E, i64 40), align 8, !dbg !301, !tbaa !302 ; /usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h:207:21
  %0 = tail call i32 @__cxa_atexit(ptr nonnull @_ZNSt3mapINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEEN6Kokkos5Tools12Experimental13TeamSizeTunerESt4lessIS5_ESaISt4pairIKS5_S9_EEED2Ev, ptr nonnull @_ZN6Kokkos5Tools12Experimental4ImplL11team_tunersB5cxx11E, ptr nonnull @__dso_handle) #12, !dbg !303 ; probe.cpp:0
  ret void
}

; Function Attrs: nocallback nofree nosync nounwind willreturn memory(inaccessiblemem: write)
declare void @llvm.assume(i1 noundef) #10

attributes #0 = { mustprogress nounwind uwtable "min-legal-vector-width"="0" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="znver3" "target-features"="+adx,+aes,+avx,+avx2,+bmi,+bmi2,+clflushopt,+clwb,+clzero,+crc32,+cx16,+cx8,+f16c,+fma,+fsgsbase,+fxsr,+invpcid,+lzcnt,+mmx,+movbe,+mwaitx,+pclmul,+pku,+popcnt,+prfchw,+rdpid,+rdpru,+rdrnd,+rdseed,+sahf,+sha,+sse,+sse2,+sse3,+sse4.1,+sse4.2,+sse4a,+ssse3,+vaes,+vpclmulqdq,+wbnoinvd,+x87,+xsave,+xsavec,+xsaveopt,+xsaves" }
attributes #1 = { nofree nounwind }
attributes #2 = { mustprogress nofree noinline nounwind willreturn memory(write) uwtable "min-legal-vector-width"="0" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="znver3" "target-features"="+adx,+aes,+avx,+avx2,+bmi,+bmi2,+clflushopt,+clwb,+clzero,+crc32,+cx16,+cx8,+f16c,+fma,+fsgsbase,+fxsr,+invpcid,+lzcnt,+mmx,+movbe,+mwaitx,+pclmul,+pku,+popcnt,+prfchw,+rdpid,+rdpru,+rdrnd,+rdseed,+sahf,+sha,+sse,+sse2,+sse3,+sse4.1,+sse4.2,+sse4a,+ssse3,+vaes,+vpclmulqdq,+wbnoinvd,+x87,+xsave,+xsavec,+xsaveopt,+xsaves" }
attributes #3 = { noinline noreturn nounwind uwtable "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="znver3" "target-features"="+adx,+aes,+avx,+avx2,+bmi,+bmi2,+clflushopt,+clwb,+clzero,+crc32,+cx16,+cx8,+f16c,+fma,+fsgsbase,+fxsr,+invpcid,+lzcnt,+mmx,+movbe,+mwaitx,+pclmul,+pku,+popcnt,+prfchw,+rdpid,+rdpru,+rdrnd,+rdseed,+sahf,+sha,+sse,+sse2,+sse3,+sse4.1,+sse4.2,+sse4a,+ssse3,+vaes,+vpclmulqdq,+wbnoinvd,+x87,+xsave,+xsavec,+xsaveopt,+xsaves" }
attributes #4 = { cold nofree noreturn }
attributes #5 = { mustprogress uwtable "min-legal-vector-width"="0" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="znver3" "target-features"="+adx,+aes,+avx,+avx2,+bmi,+bmi2,+clflushopt,+clwb,+clzero,+crc32,+cx16,+cx8,+f16c,+fma,+fsgsbase,+fxsr,+invpcid,+lzcnt,+mmx,+movbe,+mwaitx,+pclmul,+pku,+popcnt,+prfchw,+rdpid,+rdpru,+rdrnd,+rdseed,+sahf,+sha,+sse,+sse2,+sse3,+sse4.1,+sse4.2,+sse4a,+ssse3,+vaes,+vpclmulqdq,+wbnoinvd,+x87,+xsave,+xsavec,+xsaveopt,+xsaves" }
attributes #6 = { nobuiltin nounwind "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="znver3" "target-features"="+adx,+aes,+avx,+avx2,+bmi,+bmi2,+clflushopt,+clwb,+clzero,+crc32,+cx16,+cx8,+f16c,+fma,+fsgsbase,+fxsr,+invpcid,+lzcnt,+mmx,+movbe,+mwaitx,+pclmul,+pku,+popcnt,+prfchw,+rdpid,+rdpru,+rdrnd,+rdseed,+sahf,+sha,+sse,+sse2,+sse3,+sse4.1,+sse4.2,+sse4a,+ssse3,+vaes,+vpclmulqdq,+wbnoinvd,+x87,+xsave,+xsavec,+xsaveopt,+xsaves" }
attributes #7 = { mustprogress nofree nounwind willreturn memory(write) "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="znver3" "target-features"="+adx,+aes,+avx,+avx2,+bmi,+bmi2,+clflushopt,+clwb,+clzero,+crc32,+cx16,+cx8,+f16c,+fma,+fsgsbase,+fxsr,+invpcid,+lzcnt,+mmx,+movbe,+mwaitx,+pclmul,+pku,+popcnt,+prfchw,+rdpid,+rdpru,+rdrnd,+rdseed,+sahf,+sha,+sse,+sse2,+sse3,+sse4.1,+sse4.2,+sse4a,+ssse3,+vaes,+vpclmulqdq,+wbnoinvd,+x87,+xsave,+xsavec,+xsaveopt,+xsaves" }
attributes #8 = { mustprogress nocallback nofree nosync nounwind speculatable willreturn memory(none) }
attributes #9 = { nofree nounwind uwtable "min-legal-vector-width"="0" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="znver3" "target-features"="+adx,+aes,+avx,+avx2,+bmi,+bmi2,+clflushopt,+clwb,+clzero,+crc32,+cx16,+cx8,+f16c,+fma,+fsgsbase,+fxsr,+invpcid,+lzcnt,+mmx,+movbe,+mwaitx,+pclmul,+pku,+popcnt,+prfchw,+rdpid,+rdpru,+rdrnd,+rdseed,+sahf,+sha,+sse,+sse2,+sse3,+sse4.1,+sse4.2,+sse4a,+ssse3,+vaes,+vpclmulqdq,+wbnoinvd,+x87,+xsave,+xsavec,+xsaveopt,+xsaves" }
attributes #10 = { nocallback nofree nosync nounwind willreturn memory(inaccessiblemem: write) }
attributes #11 = { noreturn nounwind }
attributes #12 = { nounwind }
attributes #13 = { builtin nounwind }

!llvm.dbg.cu = !{!0}
!llvm.module.flags = !{!2, !3, !4}
!llvm.ident = !{!5}

!0 = distinct !DICompileUnit(language: DW_LANG_C_plus_plus_14, file: !1, producer: "Cray clang version 20.0.0  (95f64a80f3ee7b4a305f3014c51d36f78b24d91e)", isOptimized: true, runtimeVersion: 0, emissionKind: NoDebug, splitDebugInlining: false, nameTableKind: None)
!1 = !DIFile(filename: "/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/passive-cast-proof/probe.cpp", directory: "/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/passive-cast-proof/cpu")
!2 = !{i32 2, !"Debug Info Version", i32 3}
!3 = !{i32 1, !"wchar_size", i32 4}
!4 = !{i32 7, !"uwtable", i32 2}
!5 = !{!"Cray clang version 20.0.0  (95f64a80f3ee7b4a305f3014c51d36f78b24d91e)"}
!6 = distinct !DISubprogram(name: "~map", scope: !7, file: !7, line: 314, type: !8, scopeLine: 314, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!7 = !DIFile(filename: "/usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_map.h", directory: "")
!8 = !DISubroutineType(types: !9)
!9 = !{}
!10 = !DILocation(line: 733, column: 64, scope: !11, inlinedAt: !13)
!11 = distinct !DISubprogram(name: "_M_mbegin", scope: !12, file: !12, line: 732, type: !8, scopeLine: 733, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!12 = !DIFile(filename: "/usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_tree.h", directory: "")
!13 = distinct !DILocation(line: 737, column: 16, scope: !14, inlinedAt: !15)
!14 = distinct !DISubprogram(name: "_M_begin", scope: !12, file: !12, line: 736, type: !8, scopeLine: 737, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!15 = distinct !DILocation(line: 982, column: 18, scope: !16, inlinedAt: !17)
!16 = distinct !DISubprogram(name: "~_Rb_tree", scope: !12, file: !12, line: 981, type: !8, scopeLine: 982, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!17 = distinct !DILocation(line: 314, column: 22, scope: !6)
!18 = !{!19, !24, i64 8}
!19 = !{!"_ZTSSt15_Rb_tree_header", !20, i64 0, !26, i64 32}
!20 = !{!"_ZTSSt18_Rb_tree_node_base", !21, i64 0, !24, i64 8, !24, i64 16, !24, i64 24}
!21 = !{!"_ZTSSt14_Rb_tree_color", !22, i64 0}
!22 = !{!"omnipotent char", !23, i64 0}
!23 = !{!"Simple C++ TBAA"}
!24 = !{!"p1 _ZTSSt18_Rb_tree_node_base", !25, i64 0}
!25 = !{!"any pointer", !22, i64 0}
!26 = !{!"long", !22, i64 0}
!27 = !DILocation(line: 982, column: 9, scope: !16, inlinedAt: !17)
!28 = !DILocation(line: 314, column: 22, scope: !6)
!29 = distinct !DISubprogram(name: "wo2_passive_specific", scope: !30, file: !30, line: 11, type: !8, scopeLine: 11, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!30 = !DIFile(filename: "probe.cpp", directory: "/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/passive-cast-proof")
!31 = !DILocation(line: 18, column: 21, scope: !32, inlinedAt: !34)
!32 = distinct !DISubprogram(name: "PassiveSpecificInvariants", scope: !33, file: !33, line: 16, type: !8, scopeLine: 17, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!33 = !DIFile(filename: "include/eos/cgl_passive.hpp", directory: "/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research/passive-cast-proof/cpu")
!34 = distinct !DILocation(line: 12, column: 18, scope: !29)
!35 = !{!36, !36, i64 0}
!36 = !{!"int", !22, i64 0}
!37 = !DILocation(line: 18, column: 38, scope: !32, inlinedAt: !34)
!38 = !DILocation(line: 19, column: 11, scope: !32, inlinedAt: !34)
!39 = !DILocation(line: 19, column: 21, scope: !32, inlinedAt: !34)
!40 = !DILocation(line: 19, column: 32, scope: !32, inlinedAt: !34)
!41 = !DILocation(line: 20, column: 11, scope: !32, inlinedAt: !34)
!42 = !DILocation(line: 20, column: 24, scope: !32, inlinedAt: !34)
!43 = !DILocation(line: 20, column: 22, scope: !32, inlinedAt: !34)
!44 = !DILocation(line: 20, column: 34, scope: !32, inlinedAt: !34)
!45 = !DILocation(line: 20, column: 45, scope: !32, inlinedAt: !34)
!46 = !DILocation(line: 13, column: 10, scope: !29)
!47 = !{!48, !48, i64 0}
!48 = !{!"double", !22, i64 0}
!49 = !DILocation(line: 14, column: 3, scope: !29)
!50 = !DILocation(line: 14, column: 10, scope: !29)
!51 = !DILocation(line: 15, column: 1, scope: !29)
!52 = distinct !DISubprogram(name: "wo2_passive_encode", scope: !30, file: !30, line: 17, type: !8, scopeLine: 17, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!53 = !DILocation(line: 18, column: 21, scope: !32, inlinedAt: !54)
!54 = distinct !DILocation(line: 26, column: 25, scope: !55, inlinedAt: !56)
!55 = distinct !DISubprogram(name: "PassiveEncode", scope: !33, file: !33, line: 24, type: !8, scopeLine: 25, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!56 = distinct !DILocation(line: 18, column: 18, scope: !52)
!57 = !DILocation(line: 18, column: 38, scope: !32, inlinedAt: !54)
!58 = !DILocation(line: 19, column: 11, scope: !32, inlinedAt: !54)
!59 = !DILocation(line: 19, column: 21, scope: !32, inlinedAt: !54)
!60 = !DILocation(line: 19, column: 32, scope: !32, inlinedAt: !54)
!61 = !DILocation(line: 20, column: 11, scope: !32, inlinedAt: !54)
!62 = !DILocation(line: 20, column: 24, scope: !32, inlinedAt: !54)
!63 = !DILocation(line: 20, column: 22, scope: !32, inlinedAt: !54)
!64 = !DILocation(line: 20, column: 34, scope: !32, inlinedAt: !54)
!65 = !DILocation(line: 20, column: 45, scope: !32, inlinedAt: !54)
!66 = !DILocation(line: 27, column: 14, scope: !55, inlinedAt: !56)
!67 = !DILocation(line: 27, column: 23, scope: !55, inlinedAt: !56)
!68 = !DILocation(line: 19, column: 10, scope: !52)
!69 = !DILocation(line: 20, column: 3, scope: !52)
!70 = !DILocation(line: 20, column: 10, scope: !52)
!71 = !DILocation(line: 21, column: 1, scope: !52)
!72 = distinct !DISubprogram(name: "_M_erase", scope: !12, file: !12, line: 1927, type: !8, scopeLine: 1928, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!73 = !DILocation(line: 1930, column: 18, scope: !72)
!74 = !{i32 1}
!75 = !DILocation(line: 1930, column: 7, scope: !72)
!76 = !DILocation(line: 786, column: 45, scope: !77, inlinedAt: !78)
!77 = distinct !DISubprogram(name: "_S_right", scope: !12, file: !12, line: 785, type: !8, scopeLine: 786, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!78 = distinct !DILocation(line: 1932, column: 13, scope: !72)
!79 = !{!20, !24, i64 24}
!80 = !DILocation(line: 1932, column: 4, scope: !72)
!81 = !DILocation(line: 778, column: 45, scope: !82, inlinedAt: !83)
!82 = distinct !DISubprogram(name: "_S_left", scope: !12, file: !12, line: 777, type: !8, scopeLine: 778, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!83 = distinct !DILocation(line: 1933, column: 21, scope: !72)
!84 = !{!20, !24, i64 16}
!85 = !DILocation(line: 231, column: 16, scope: !86, inlinedAt: !87)
!86 = distinct !DISubprogram(name: "_M_valptr", scope: !12, file: !12, line: 230, type: !8, scopeLine: 231, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!87 = distinct !DILocation(line: 621, column: 55, scope: !88, inlinedAt: !89)
!88 = distinct !DISubprogram(name: "_M_destroy_node", scope: !12, file: !12, line: 616, type: !8, scopeLine: 617, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!89 = distinct !DILocation(line: 629, column: 2, scope: !90, inlinedAt: !91)
!90 = distinct !DISubprogram(name: "_M_drop_node", scope: !12, file: !12, line: 627, type: !8, scopeLine: 628, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!91 = distinct !DILocation(line: 1934, column: 4, scope: !72)
!92 = !DILocation(line: 279, column: 7, scope: !93, inlinedAt: !95)
!93 = distinct !DISubprogram(name: "~MultidimensionalSparseTuningProblem", scope: !94, file: !94, line: 279, type: !8, scopeLine: 279, flags: DIFlagArtificial | DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!94 = !DIFile(filename: "source-unfused/kokkos/core/src/Kokkos_Tuners.hpp", directory: "/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research")
!95 = distinct !DILocation(line: 406, column: 7, scope: !96, inlinedAt: !97)
!96 = distinct !DISubprogram(name: "~TeamSizeTuner", scope: !94, file: !94, line: 406, type: !8, scopeLine: 406, flags: DIFlagArtificial | DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!97 = distinct !DILocation(line: 284, column: 12, scope: !98, inlinedAt: !100)
!98 = distinct !DISubprogram(name: "~pair", scope: !99, file: !99, line: 284, type: !8, scopeLine: 284, flags: DIFlagArtificial | DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!99 = !DIFile(filename: "/usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_pair.h", directory: "")
!100 = distinct !DILocation(line: 198, column: 10, scope: !101, inlinedAt: !103)
!101 = distinct !DISubprogram(name: "destroy<std::pair<const std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Kokkos::Tools::Experimental::TeamSizeTuner> >", scope: !102, file: !102, line: 196, type: !8, scopeLine: 198, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!102 = !DIFile(filename: "/usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/new_allocator.h", directory: "")
!103 = distinct !DILocation(line: 554, column: 8, scope: !104, inlinedAt: !106)
!104 = distinct !DISubprogram(name: "destroy<std::pair<const std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, Kokkos::Tools::Experimental::TeamSizeTuner> >", scope: !105, file: !105, line: 550, type: !8, scopeLine: 552, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!105 = !DIFile(filename: "/usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/alloc_traits.h", directory: "")
!106 = distinct !DILocation(line: 621, column: 2, scope: !88, inlinedAt: !89)
!107 = !DILocation(line: 735, column: 30, scope: !108, inlinedAt: !110)
!108 = distinct !DISubprogram(name: "~vector", scope: !109, file: !109, line: 733, type: !8, scopeLine: 734, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!109 = !DIFile(filename: "/usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_vector.h", directory: "")
!110 = distinct !DILocation(line: 279, column: 7, scope: !93, inlinedAt: !95)
!111 = !{!112, !113, i64 0}
!112 = !{!"_ZTSNSt12_Vector_baseINSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEESaIS5_EE17_Vector_impl_dataE", !113, i64 0, !113, i64 8, !113, i64 16}
!113 = !{!"p1 _ZTSNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE", !25, i64 0}
!114 = !DILocation(line: 735, column: 54, scope: !108, inlinedAt: !110)
!115 = !{!112, !113, i64 8}
!116 = !DILocation(line: 162, column: 19, scope: !117, inlinedAt: !119)
!117 = distinct !DISubprogram(name: "__destroy<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > *>", scope: !118, file: !118, line: 160, type: !8, scopeLine: 161, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!118 = !DIFile(filename: "/usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/stl_construct.h", directory: "")
!119 = distinct !DILocation(line: 195, column: 7, scope: !120, inlinedAt: !121)
!120 = distinct !DISubprogram(name: "_Destroy<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > *>", scope: !118, file: !118, line: 182, type: !8, scopeLine: 183, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!121 = distinct !DILocation(line: 944, column: 7, scope: !122, inlinedAt: !123)
!122 = distinct !DISubprogram(name: "_Destroy<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > *, std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >", scope: !105, file: !105, line: 941, type: !8, scopeLine: 943, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!123 = distinct !DILocation(line: 735, column: 2, scope: !108, inlinedAt: !110)
!124 = !DILocation(line: 162, column: 4, scope: !117, inlinedAt: !119)
!125 = !DILocation(line: 228, column: 28, scope: !126, inlinedAt: !128)
!126 = distinct !DISubprogram(name: "_M_data", scope: !127, file: !127, line: 227, type: !8, scopeLine: 228, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!127 = !DIFile(filename: "/usr/lib64/gcc/x86_64-suse-linux/14/../../../../include/c++/14/bits/basic_string.h", directory: "")
!128 = distinct !DILocation(line: 269, column: 6, scope: !129, inlinedAt: !130)
!129 = distinct !DISubprogram(name: "_M_is_local", scope: !127, file: !127, line: 267, type: !8, scopeLine: 268, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!130 = distinct !DILocation(line: 287, column: 7, scope: !131, inlinedAt: !132)
!131 = distinct !DISubprogram(name: "_M_dispose", scope: !127, file: !127, line: 285, type: !8, scopeLine: 286, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!132 = distinct !DILocation(line: 809, column: 9, scope: !133, inlinedAt: !134)
!133 = distinct !DISubprogram(name: "~basic_string", scope: !127, file: !127, line: 808, type: !8, scopeLine: 809, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!134 = distinct !DILocation(line: 151, column: 19, scope: !135, inlinedAt: !136)
!135 = distinct !DISubprogram(name: "_Destroy<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > >", scope: !118, file: !118, line: 146, type: !8, scopeLine: 147, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!136 = distinct !DILocation(line: 163, column: 6, scope: !117, inlinedAt: !119)
!137 = !{!138, !140, i64 0}
!138 = !{!"_ZTSNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE", !139, i64 0, !26, i64 8, !22, i64 16}
!139 = !{!"_ZTSNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEE12_Alloc_hiderE", !140, i64 0}
!140 = !{!"p1 omnipotent char", !25, i64 0}
!141 = !DILocation(line: 246, column: 57, scope: !142, inlinedAt: !143)
!142 = distinct !DISubprogram(name: "_M_local_data", scope: !127, file: !127, line: 243, type: !8, scopeLine: 244, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!143 = distinct !DILocation(line: 269, column: 19, scope: !129, inlinedAt: !130)
!144 = !DILocation(line: 269, column: 16, scope: !129, inlinedAt: !130)
!145 = !DILocation(line: 271, column: 10, scope: !129, inlinedAt: !130)
!146 = !{!138, !26, i64 8}
!147 = !DILocation(line: 271, column: 27, scope: !129, inlinedAt: !130)
!148 = !DILocation(line: 287, column: 6, scope: !131, inlinedAt: !132)
!149 = !DILocation(line: 288, column: 15, scope: !131, inlinedAt: !132)
!150 = !{!22, !22, i64 0}
!151 = !DILocation(line: 294, column: 73, scope: !152, inlinedAt: !153)
!152 = distinct !DISubprogram(name: "_M_destroy", scope: !127, file: !127, line: 293, type: !8, scopeLine: 294, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!153 = distinct !DILocation(line: 288, column: 4, scope: !131, inlinedAt: !132)
!154 = !DILocation(line: 172, column: 2, scope: !155, inlinedAt: !156)
!155 = distinct !DISubprogram(name: "deallocate", scope: !102, file: !102, line: 156, type: !8, scopeLine: 157, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!156 = distinct !DILocation(line: 513, column: 13, scope: !157, inlinedAt: !158)
!157 = distinct !DISubprogram(name: "deallocate", scope: !105, file: !105, line: 512, type: !8, scopeLine: 513, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!158 = distinct !DILocation(line: 294, column: 9, scope: !152, inlinedAt: !153)
!159 = !DILocation(line: 288, column: 4, scope: !131, inlinedAt: !132)
!160 = !DILocation(line: 162, column: 30, scope: !117, inlinedAt: !119)
!161 = distinct !{!161, !124, !162, !163}
!162 = !DILocation(line: 163, column: 46, scope: !117, inlinedAt: !119)
!163 = !{!"llvm.loop.mustprogress"}
!164 = !DILocation(line: 368, column: 24, scope: !165, inlinedAt: !166)
!165 = distinct !DISubprogram(name: "~_Vector_base", scope: !109, file: !109, line: 366, type: !8, scopeLine: 367, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!166 = distinct !DILocation(line: 738, column: 7, scope: !108, inlinedAt: !110)
!167 = !DILocation(line: 388, column: 6, scope: !168, inlinedAt: !169)
!168 = distinct !DISubprogram(name: "_M_deallocate", scope: !109, file: !109, line: 385, type: !8, scopeLine: 386, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!169 = distinct !DILocation(line: 368, column: 2, scope: !165, inlinedAt: !166)
!170 = !DILocation(line: 369, column: 17, scope: !165, inlinedAt: !166)
!171 = !{!112, !113, i64 16}
!172 = !DILocation(line: 369, column: 35, scope: !165, inlinedAt: !166)
!173 = !DILocation(line: 172, column: 2, scope: !174, inlinedAt: !175)
!174 = distinct !DISubprogram(name: "deallocate", scope: !102, file: !102, line: 156, type: !8, scopeLine: 157, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!175 = distinct !DILocation(line: 513, column: 13, scope: !176, inlinedAt: !177)
!176 = distinct !DISubprogram(name: "deallocate", scope: !105, file: !105, line: 512, type: !8, scopeLine: 513, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!177 = distinct !DILocation(line: 389, column: 4, scope: !168, inlinedAt: !169)
!178 = !DILocation(line: 389, column: 4, scope: !168, inlinedAt: !169)
!179 = !DILocation(line: 69, column: 8, scope: !180, inlinedAt: !181)
!180 = distinct !DISubprogram(name: "~ValueHierarchyNode", scope: !94, file: !94, line: 69, type: !8, scopeLine: 69, flags: DIFlagArtificial | DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!181 = distinct !DILocation(line: 279, column: 7, scope: !93, inlinedAt: !95)
!182 = !DILocation(line: 735, column: 30, scope: !183, inlinedAt: !184)
!183 = distinct !DISubprogram(name: "~vector", scope: !109, file: !109, line: 733, type: !8, scopeLine: 734, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!184 = distinct !DILocation(line: 69, column: 8, scope: !180, inlinedAt: !181)
!185 = !{!186, !187, i64 0}
!186 = !{!"_ZTSNSt12_Vector_baseIN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEESaIS5_EE17_Vector_impl_dataE", !187, i64 0, !187, i64 8, !187, i64 16}
!187 = !{!"p1 _ZTSN6Kokkos5Tools12Experimental4Impl18ValueHierarchyNodeIlvEE", !25, i64 0}
!188 = !DILocation(line: 735, column: 54, scope: !183, inlinedAt: !184)
!189 = !{!186, !187, i64 8}
!190 = !DILocation(line: 162, column: 19, scope: !191, inlinedAt: !192)
!191 = distinct !DISubprogram(name: "__destroy<Kokkos::Tools::Experimental::Impl::ValueHierarchyNode<long, void> *>", scope: !118, file: !118, line: 160, type: !8, scopeLine: 161, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!192 = distinct !DILocation(line: 195, column: 7, scope: !193, inlinedAt: !194)
!193 = distinct !DISubprogram(name: "_Destroy<Kokkos::Tools::Experimental::Impl::ValueHierarchyNode<long, void> *>", scope: !118, file: !118, line: 182, type: !8, scopeLine: 183, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!194 = distinct !DILocation(line: 944, column: 7, scope: !195, inlinedAt: !196)
!195 = distinct !DISubprogram(name: "_Destroy<Kokkos::Tools::Experimental::Impl::ValueHierarchyNode<long, void> *, Kokkos::Tools::Experimental::Impl::ValueHierarchyNode<long, void> >", scope: !105, file: !105, line: 941, type: !8, scopeLine: 943, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!196 = distinct !DILocation(line: 735, column: 2, scope: !183, inlinedAt: !184)
!197 = !DILocation(line: 162, column: 4, scope: !191, inlinedAt: !192)
!198 = !DILocation(line: 368, column: 24, scope: !199, inlinedAt: !200)
!199 = distinct !DISubprogram(name: "~_Vector_base", scope: !109, file: !109, line: 366, type: !8, scopeLine: 367, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!200 = distinct !DILocation(line: 738, column: 7, scope: !201, inlinedAt: !202)
!201 = distinct !DISubprogram(name: "~vector", scope: !109, file: !109, line: 733, type: !8, scopeLine: 734, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!202 = distinct !DILocation(line: 69, column: 8, scope: !203, inlinedAt: !204)
!203 = distinct !DISubprogram(name: "~ValueHierarchyNode", scope: !94, file: !94, line: 69, type: !8, scopeLine: 69, flags: DIFlagArtificial | DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!204 = distinct !DILocation(line: 151, column: 19, scope: !205, inlinedAt: !206)
!205 = distinct !DISubprogram(name: "_Destroy<Kokkos::Tools::Experimental::Impl::ValueHierarchyNode<long, void> >", scope: !118, file: !118, line: 146, type: !8, scopeLine: 147, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!206 = distinct !DILocation(line: 163, column: 6, scope: !191, inlinedAt: !192)
!207 = !{!208, !209, i64 0}
!208 = !{!"_ZTSNSt12_Vector_baseIlSaIlEE17_Vector_impl_dataE", !209, i64 0, !209, i64 8, !209, i64 16}
!209 = !{!"p1 long", !25, i64 0}
!210 = !DILocation(line: 388, column: 6, scope: !211, inlinedAt: !212)
!211 = distinct !DISubprogram(name: "_M_deallocate", scope: !109, file: !109, line: 385, type: !8, scopeLine: 386, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!212 = distinct !DILocation(line: 368, column: 2, scope: !199, inlinedAt: !200)
!213 = !DILocation(line: 369, column: 17, scope: !199, inlinedAt: !200)
!214 = !{!208, !209, i64 16}
!215 = !DILocation(line: 369, column: 35, scope: !199, inlinedAt: !200)
!216 = !DILocation(line: 172, column: 2, scope: !217, inlinedAt: !218)
!217 = distinct !DISubprogram(name: "deallocate", scope: !102, file: !102, line: 156, type: !8, scopeLine: 157, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!218 = distinct !DILocation(line: 513, column: 13, scope: !219, inlinedAt: !220)
!219 = distinct !DISubprogram(name: "deallocate", scope: !105, file: !105, line: 512, type: !8, scopeLine: 513, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!220 = distinct !DILocation(line: 389, column: 4, scope: !211, inlinedAt: !212)
!221 = !DILocation(line: 389, column: 4, scope: !211, inlinedAt: !212)
!222 = !DILocation(line: 162, column: 30, scope: !191, inlinedAt: !192)
!223 = distinct !{!223, !197, !224, !163}
!224 = !DILocation(line: 163, column: 46, scope: !191, inlinedAt: !192)
!225 = !DILocation(line: 368, column: 24, scope: !226, inlinedAt: !227)
!226 = distinct !DISubprogram(name: "~_Vector_base", scope: !109, file: !109, line: 366, type: !8, scopeLine: 367, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!227 = distinct !DILocation(line: 738, column: 7, scope: !183, inlinedAt: !184)
!228 = !DILocation(line: 388, column: 6, scope: !229, inlinedAt: !230)
!229 = distinct !DISubprogram(name: "_M_deallocate", scope: !109, file: !109, line: 385, type: !8, scopeLine: 386, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!230 = distinct !DILocation(line: 368, column: 2, scope: !226, inlinedAt: !227)
!231 = !DILocation(line: 369, column: 17, scope: !226, inlinedAt: !227)
!232 = !{!186, !187, i64 16}
!233 = !DILocation(line: 369, column: 35, scope: !226, inlinedAt: !227)
!234 = !DILocation(line: 172, column: 2, scope: !235, inlinedAt: !236)
!235 = distinct !DISubprogram(name: "deallocate", scope: !102, file: !102, line: 156, type: !8, scopeLine: 157, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!236 = distinct !DILocation(line: 513, column: 13, scope: !237, inlinedAt: !238)
!237 = distinct !DISubprogram(name: "deallocate", scope: !105, file: !105, line: 512, type: !8, scopeLine: 513, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!238 = distinct !DILocation(line: 389, column: 4, scope: !229, inlinedAt: !230)
!239 = !DILocation(line: 389, column: 4, scope: !229, inlinedAt: !230)
!240 = !DILocation(line: 368, column: 24, scope: !199, inlinedAt: !241)
!241 = distinct !DILocation(line: 738, column: 7, scope: !201, inlinedAt: !242)
!242 = distinct !DILocation(line: 69, column: 8, scope: !180, inlinedAt: !181)
!243 = !DILocation(line: 388, column: 6, scope: !211, inlinedAt: !244)
!244 = distinct !DILocation(line: 368, column: 2, scope: !199, inlinedAt: !241)
!245 = !DILocation(line: 369, column: 17, scope: !199, inlinedAt: !241)
!246 = !DILocation(line: 369, column: 35, scope: !199, inlinedAt: !241)
!247 = !DILocation(line: 172, column: 2, scope: !217, inlinedAt: !248)
!248 = distinct !DILocation(line: 513, column: 13, scope: !219, inlinedAt: !249)
!249 = distinct !DILocation(line: 389, column: 4, scope: !211, inlinedAt: !244)
!250 = !DILocation(line: 389, column: 4, scope: !211, inlinedAt: !244)
!251 = !DILocation(line: 228, column: 28, scope: !126, inlinedAt: !252)
!252 = distinct !DILocation(line: 269, column: 6, scope: !129, inlinedAt: !253)
!253 = distinct !DILocation(line: 287, column: 7, scope: !131, inlinedAt: !254)
!254 = distinct !DILocation(line: 809, column: 9, scope: !133, inlinedAt: !255)
!255 = distinct !DILocation(line: 284, column: 12, scope: !98, inlinedAt: !100)
!256 = !DILocation(line: 246, column: 57, scope: !142, inlinedAt: !257)
!257 = distinct !DILocation(line: 269, column: 19, scope: !129, inlinedAt: !253)
!258 = !DILocation(line: 269, column: 16, scope: !129, inlinedAt: !253)
!259 = !DILocation(line: 271, column: 10, scope: !129, inlinedAt: !253)
!260 = !DILocation(line: 271, column: 27, scope: !129, inlinedAt: !253)
!261 = !DILocation(line: 287, column: 6, scope: !131, inlinedAt: !254)
!262 = !DILocation(line: 288, column: 15, scope: !131, inlinedAt: !254)
!263 = !DILocation(line: 294, column: 73, scope: !152, inlinedAt: !264)
!264 = distinct !DILocation(line: 288, column: 4, scope: !131, inlinedAt: !254)
!265 = !DILocation(line: 172, column: 2, scope: !155, inlinedAt: !266)
!266 = distinct !DILocation(line: 513, column: 13, scope: !157, inlinedAt: !267)
!267 = distinct !DILocation(line: 294, column: 9, scope: !152, inlinedAt: !264)
!268 = !DILocation(line: 288, column: 4, scope: !131, inlinedAt: !254)
!269 = !DILocation(line: 172, column: 2, scope: !270, inlinedAt: !271)
!270 = distinct !DISubprogram(name: "deallocate", scope: !102, file: !102, line: 156, type: !8, scopeLine: 157, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!271 = distinct !DILocation(line: 513, column: 13, scope: !272, inlinedAt: !273)
!272 = distinct !DISubprogram(name: "deallocate", scope: !105, file: !105, line: 512, type: !8, scopeLine: 513, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!273 = distinct !DILocation(line: 563, column: 9, scope: !274, inlinedAt: !275)
!274 = distinct !DISubprogram(name: "_M_put_node", scope: !12, file: !12, line: 562, type: !8, scopeLine: 563, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!275 = distinct !DILocation(line: 630, column: 2, scope: !90, inlinedAt: !91)
!276 = distinct !{!276, !75, !277, !163}
!277 = !DILocation(line: 1936, column: 2, scope: !72)
!278 = !DILocation(line: 1937, column: 5, scope: !72)
!279 = distinct !DISubprogram(linkageName: "_GLOBAL__sub_I_probe.cpp", scope: !30, file: !30, type: !8, flags: DIFlagArtificial, spFlags: DISPFlagLocalToUnit | DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!280 = !DILocation(line: 171, column: 26, scope: !281, inlinedAt: !282)
!281 = distinct !DISubprogram(name: "_Rb_tree_header", scope: !12, file: !12, line: 169, type: !8, scopeLine: 170, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!282 = distinct !DILocation(line: 665, column: 4, scope: !283, inlinedAt: !284)
!283 = distinct !DISubprogram(name: "_Rb_tree_impl", scope: !12, file: !12, line: 665, type: !8, scopeLine: 670, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!284 = distinct !DILocation(line: 926, column: 7, scope: !285, inlinedAt: !286)
!285 = distinct !DISubprogram(name: "_Rb_tree", scope: !12, file: !12, line: 926, type: !8, scopeLine: 926, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!286 = distinct !DILocation(line: 197, column: 7, scope: !287, inlinedAt: !288)
!287 = distinct !DISubprogram(name: "map", scope: !7, file: !7, line: 197, type: !8, scopeLine: 197, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!288 = distinct !DILocation(line: 37, column: 5, scope: !289, inlinedAt: !292)
!289 = !DILexicalBlockFile(scope: !291, file: !290, discriminator: 0)
!290 = !DIFile(filename: "source-unfused/kokkos/core/src/impl/Kokkos_Tools_Generic.hpp", directory: "/lustre/orion/ast207/proj-shared/dfielding/CGL/WO2/final-research")
!291 = distinct !DISubprogram(name: "__cxx_global_var_init", scope: !30, file: !30, type: !8, flags: DIFlagArtificial, spFlags: DISPFlagLocalToUnit | DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!292 = distinct !DILocation(line: 0, scope: !279)
!293 = !{!19, !21, i64 0}
!294 = !DILocation(line: 204, column: 27, scope: !295, inlinedAt: !296)
!295 = distinct !DISubprogram(name: "_M_reset", scope: !12, file: !12, line: 202, type: !8, scopeLine: 203, flags: DIFlagPrototyped, spFlags: DISPFlagDefinition | DISPFlagOptimized, unit: !0)
!296 = distinct !DILocation(line: 172, column: 7, scope: !281, inlinedAt: !282)
!297 = !DILocation(line: 205, column: 25, scope: !295, inlinedAt: !296)
!298 = !{!19, !24, i64 16}
!299 = !DILocation(line: 206, column: 26, scope: !295, inlinedAt: !296)
!300 = !{!19, !24, i64 24}
!301 = !DILocation(line: 207, column: 21, scope: !295, inlinedAt: !296)
!302 = !{!19, !26, i64 32}
!303 = !DILocation(line: 0, scope: !291, inlinedAt: !292)
