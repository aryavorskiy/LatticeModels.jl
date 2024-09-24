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
    @show rep
    @test_broken length(JET.get_reports(rep)) == 0
end
