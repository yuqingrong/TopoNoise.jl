using Test,TopoNoise
@testset "H+5 direct-projector scout and production CLI" begin
    path=joinpath(@__DIR__,"..","..","examples","compare_toric3d_vertex_projector_feedback.jl")
    @test isfile(path)
    if isfile(path)
        mod=Module(:H5Example);Base.include(mod,path)
        parser=Base.invokelatest(getproperty,mod,:parse_vertex_args)
        args=Base.invokelatest(parser,["--stage","scout"])
        @test args["dimensions"]=="3x3x3,5x5x5,7x7x7"
        @test parse(Int,args["shots"])==500 && args["seed"]=="1234"
        @test parse.(Float64,split(args["rates"],','))==[0.,1e-4,3e-4,1e-3,3e-3,.01,.02,.03,.05,.075,.10]
        @test args["protocols"]=="none,one_step,window_complete,full_prefix"
        @test_throws ArgumentError Base.invokelatest(parser,["--protocols","none"])
        @test_throws ArgumentError Base.invokelatest(parser,["--stage","production"])
        @test_throws ArgumentError Base.invokelatest(parser,["--bogus","1"])
    end
end
@testset "H+5 selection schema supports separate channel grids" begin
    @test isdefined(TopoNoise,:Toric3DVertexProductionSelection)
    if isdefined(TopoNoise,:Toric3DVertexProductionSelection)
        dims=[[i,i,i] for i in (3,5,7,9)];ps=[.008,.009,.01,.011,.012]
        ch=Dict("dimensions_T_Lu_Lv"=>dims,"rates"=>ps,"fit_windows"=>Dict(string(p)=>[[.008,.012]] for p in TORIC3D_VERTEX_PROTOCOLS),"candidate_centers"=>Dict(string(p)=>[.01] for p in TORIC3D_VERTEX_PROTOCOLS),"extension_states"=>Dict(string(p)=>"resolved_candidates" for p in TORIC3D_VERTEX_PROTOCOLS),"flags"=>Dict(string(p)=>Dict("threshold_ineligible"=>false) for p in TORIC3D_VERTEX_PROTOCOLS),"saturation_bounds"=>[1e-8,.49],"timing_model"=>Dict("shot"=>[0.,1.],"build"=>[0.,1.]),"projected_seconds"=>100.)
        payload=Dict("format_version"=>1,"rule_version"=>"h5-scout-selection-v1","status"=>"ready","extension_round"=>0,"scout_inputs"=>[Dict("csv_sha256"=>repeat("a",64),"metadata_sha256"=>repeat("b",64),"path_ids"=>["fixture.csv","fixture.toml"],"stage"=>"scout","extension_round"=>0,"parent_selection_sha256"=>"")],"seed"=>1234,"decoder_settings"=>TopoNoise._vertex_settings(Toric3DDecoderSettings()),"source_hashes"=>TopoNoise._vertex_source_hashes(),"production_shots"=>1000,"extension_grids"=>Dict(),"channels"=>Dict("x_only"=>deepcopy(ch),"z_only"=>deepcopy(ch)))
        payload["channels"]["z_only"]["rates"]=[.01,.02,.03,.04,.05]
        for p in TORIC3D_VERTEX_PROTOCOLS
            payload["channels"]["z_only"]["fit_windows"][string(p)]=[[.01,.05]]
            payload["channels"]["z_only"]["candidate_centers"][string(p)]=[.03]
        end
        mktempdir() do dir
            path=joinpath(dir,"selection.json");TopoNoise._vertex_write_json(path,payload)
            sel=load_toric3d_vertex_production_selection(path)
            @test sel.channel_grids[:x_only].error_rates!=sel.channel_grids[:z_only].error_rates
            high=deepcopy(payload);high["channels"]["x_only"]["rates"]=[.352,.396,.44,.484,.528]
            for protocol in string.(TORIC3D_VERTEX_PROTOCOLS)
                high["channels"]["x_only"]["candidate_centers"][protocol]=[.44]
                high["channels"]["x_only"]["fit_windows"][protocol]=[[.352,.528]]
            end
            TopoNoise._vertex_write_json(path,high)
            @test maximum(load_toric3d_vertex_production_selection(path).channel_grids[:x_only].error_rates)==.528

            for key in ("scout_inputs","decoder_settings","rule_version","source_hashes")
                bad=deepcopy(payload);delete!(bad,key);TopoNoise._vertex_write_json(path,bad)
                @test_throws ArgumentError load_toric3d_vertex_production_selection(path)
            end
            bad=deepcopy(payload);bad["channels"]["x_only"]["fit_windows"]["none"]=[];TopoNoise._vertex_write_json(path,bad)
            @test_throws ArgumentError load_toric3d_vertex_production_selection(path)
        end
    end
end
