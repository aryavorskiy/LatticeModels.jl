using JET
import LinearAlgebra, ProgressMeter, KrylovKit

@testset "JET" begin
    jet_ignored_modules=(
        AnyFrameModule(LinearAlgebra),
        AnyFrameModule(ProgressMeter),
        AnyFrameModule(KrylovKit),
    )
    rep = report_package(LatticeModels; toplevel_logger=nothing,
        ignored_modules=jet_ignored_modules)
    nreports = length(JET.get_reports(rep))
    MAX_REPORTS = 17
    (ENV["SHOW_JET_REPORT"] == "true" || nreports > MAX_REPORTS) && show_report(rep)
    @test nreports <= MAX_REPORTS
    @test_broken length(JET.get_reports(rep)) == 0
end
