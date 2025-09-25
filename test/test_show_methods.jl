# Comprehensive test suite for show methods
using SphericalSpatialTrees: SphericalSpatialTrees as SST
using Test
import DimensionalData as DD

@testset "Comprehensive Show Methods Test Suite" begin

    @testset "Integration Tests - ProjectionSource with different types and dimensions" begin

        @testset "Float64 with RegularGridTree - Standard dimensions" begin
            # Create test data
            x = range(-180, 180, length=361)
            y = range(-90, 90, length=181)
            data = rand(Float64, 360, 180)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            ps = SST.ProjectionSource(SST.RegularGridTree, ar)

            # Test output
            io = IOBuffer()
            show(io, MIME"text/plain"(), ps)
            output = String(take!(io))
            
            # 
            @test output == "ProjectionSource{Float64}(360×180, RegularGridTree(361×181))"

            # Verify conciseness - should be much shorter than default
            io_default = IOBuffer()
            show(io_default, ps)
            default_output = String(take!(io_default))
            @test length(output) < length(default_output) / 3  # At least 3x more concise
        end

        @testset "Float32 with RegularGridTree - Large dimensions" begin
            # Test with large arrays
            x = range(0, 100, length=1001)
            y = range(0, 50, length=501)
            data = rand(Float32, 1000, 500)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            ps = SST.ProjectionSource(SST.RegularGridTree, ar)

            io = IOBuffer()
            show(io, MIME"text/plain"(), ps)
            output = String(take!(io))

            expected = "ProjectionSource{Float32}(1000×500, RegularGridTree(1001×501))"
            @test output == expected
        end

        @testset "Int32 with RegularGridTree - Small dimensions" begin
            # Test with small arrays
            x = range(0, 1, length=11)
            y = range(0, 1, length=6)
            data = rand(Int32, 10, 5)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            ps = SST.ProjectionSource(SST.RegularGridTree, ar)

            io = IOBuffer()
            show(io, MIME"text/plain"(), ps)
            output = String(take!(io))

            expected = "ProjectionSource{Int32}(10×5, RegularGridTree(11×6))"
            @test output == expected
        end

        @testset "Edge case - Very small array (1x1)" begin
            x = range(0, 1, length=2)
            y = range(0, 1, length=2)
            data = rand(Float64, 1, 1)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            ps = SST.ProjectionSource(SST.RegularGridTree, ar)

            io = IOBuffer()
            show(io, MIME"text/plain"(), ps)
            output = String(take!(io))

            expected = "ProjectionSource{Float64}(1×1, RegularGridTree(2×2))"
            @test output == expected
        end
    end

    @testset "Integration Tests - LazyProjectedDiskArray with different configurations" begin

        @testset "Float64 with ISEACircleTree - Resolution 4" begin
            # Create source
            x = range(-180, 180, length=257)
            y = range(-90, 90, length=257)
            data = rand(Float64, 256, 256)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            source = SST.ProjectionSource(SST.RegularGridTree, ar)

            # Create target with ISEACircleTree
            target = SphericalSpatialTrees.ProjectionTarget(SST.ISEACircleTree, 4, 2)
            lpda = SphericalSpatialTrees.LazyProjectedDiskArray(source, target)

            io = IOBuffer()
            show(io, MIME"text/plain"(), lpda)
            output = String(take!(io))

            expected = "16×16×10 LazyProjectedDiskArray{Float64}"
            @test output == expected

            # Verify conciseness
            io_default = IOBuffer()
            show(io_default, lpda)
            default_output = String(take!(io_default))
            @test length(output) < length(default_output) / 3
        end

        @testset "Float32 with ISEACircleTree - Resolution 6" begin
            x = range(0, 10, length=101)
            y = range(0, 5, length=51)
            data = rand(Float32, 100, 50)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            source = SST.ProjectionSource(SST.RegularGridTree, ar)

            target = SphericalSpatialTrees.ProjectionTarget(SST.ISEACircleTree, 6, 3)
            lpda = SphericalSpatialTrees.LazyProjectedDiskArray(source, target)

            io = IOBuffer()
            show(io, MIME"text/plain"(), lpda)
            output = String(take!(io))

            expected = "64×64×10 LazyProjectedDiskArray{Float32}"
            @test output == expected
        end

        @testset "Int16 with ISEACircleTree - Resolution 2 (small)" begin
            x = range(-1, 1, length=21)
            y = range(-1, 1, length=21)
            data = rand(Int16, 20, 20)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            source = SST.ProjectionSource(SST.RegularGridTree, ar)

            target = SST.ProjectionTarget(SST.ISEACircleTree, 2, 1)
            lpda = SST.LazyProjectedDiskArray(source, target)

            io = IOBuffer()
            show(io, MIME"text/plain"(), lpda)
            output = String(take!(io))

            expected = "4×4×10 LazyProjectedDiskArray{Int16}"
            @test output == expected
        end

        @testset "Edge case - High resolution (Resolution 8)" begin
            x = range(-180, 180, length=1001)
            y = range(-90, 90, length=501)
            data = rand(Float64, 1000, 500)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            source = SST.ProjectionSource(SST.RegularGridTree, ar)

            target = SST.ProjectionTarget(SST.ISEACircleTree, 8, 4)
            lpda = SST.LazyProjectedDiskArray(source, target)

            io = IOBuffer()
            show(io, MIME"text/plain"(), lpda)
            output = String(take!(io))

            expected = "256×256×10 LazyProjectedDiskArray{Float64}"
            @test output == expected
        end
    end

    @testset "Tree Integration Tests" begin

        @testset "RegularGridTree show method consistency" begin
            # Test various grid sizes
            test_cases = [
                (1:10, 1:5, "RegularGridTree(10×5 array, chunksize)\ndimensions: 9×4"),
                (1:100, 1:50, "RegularGridTree(100×50 array, chunksize)\ndimensions: 99×49"),
                (range(0, 1, length=1001), range(0, 1, length=501), "RegularGridTree(1001×501 array, chunksize)\ndimensions: 1000×500"),
                (1:2, 1:2, "RegularGridTree(2×2 array, chunksize)\ndimensions: 1×1")  # Edge case
            ]

            for (x, y, expected) in test_cases
                tree = SST.RegularGridTree(x, y)
                io = IOBuffer()
                show(io, tree)
                output = String(take!(io))
                @test output == expected
            end
        end

        @testset "RegularGridTree copy-pastable constructor format" begin
            # Test that first line matches copy-pastable constructor pattern
            x = range(-180, 180, length=361)
            y = range(-90, 90, length=181)
            tree = SST.RegularGridTree(x, y)
            
            io = IOBuffer()
            show(io, tree)
            output = String(take!(io))
            
            # Split output into lines
            lines = split(output, '\n')
            
            # First line should be copy-pastable constructor format
            @test lines[1] == "RegularGridTree(361×181 array, chunksize)"
            
            # Second line should contain dimensions info
            @test lines[2] == "dimensions: 360×180"
            
            # Test that first line follows the expected pattern
            @test occursin(r"RegularGridTree\(\d+×\d+ array, chunksize\)", lines[1])
        end

        @testset "ISEACircleTree show method consistency" begin
            # Test various resolutions with new format
            test_cases = [
                (0, "10-leaf ISEACircleTree(0)"),
                (2, "160-leaf ISEACircleTree(2)"),
                (4, "2560-leaf ISEACircleTree(4)"),
                (6, "40960-leaf ISEACircleTree(6)"),
                (8, "655360-leaf ISEACircleTree(8)")
            ]

            for (resolution, expected) in test_cases
                tree = SST.ISEACircleTree(resolution)
                io = IOBuffer()
                show(io, tree)
                output = String(take!(io))
                @test output == expected
            end
        end

        @testset "ISEACircleTree copy-pastable constructor format" begin
            # Test that output matches the new format pattern
            tree = SST.ISEACircleTree(4)
            
            io = IOBuffer()
            show(io, tree)
            output = String(take!(io))
            
            # Should be single line with format "$(nelem)-leaf ISEACircleTree(n)"
            @test output == "2560-leaf ISEACircleTree(4)"
            
            # Test that output follows the expected pattern
            @test occursin(r"\d+-leaf ISEACircleTree\(\d+\)", output)
            
            # Test different resolutions
            for resolution in [0, 1, 2, 5, 7]
                tree = SST.ISEACircleTree(resolution)
                io = IOBuffer()
                show(io, tree)
                output = String(take!(io))
                
                # Should match the pattern with correct element count
                expected_elements = 10 * 4^resolution
                expected_output = "$(expected_elements)-leaf ISEACircleTree($resolution)"
                @test output == expected_output
                
                # Verify pattern match
                @test occursin(r"\d+-leaf ISEACircleTree\(\d+\)", output)
            end
        end

        @testset "Tree integration in ProjectionSource" begin
            # Test that tree display integrates properly
            x = range(0, 10, length=51)
            y = range(0, 5, length=26)
            data = rand(Float64, 50, 25)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            ps = SST.ProjectionSource(SST.RegularGridTree, ar)

            io = IOBuffer()
            show(io, MIME"text/plain"(), ps)
            output = String(take!(io))

            # Should contain both the ProjectionSource info and the tree info
            @test contains(output, "ProjectionSource{Float64}")
            @test contains(output, "50×25")
            @test contains(output, "RegularGridTree(51×26)")

            expected = "ProjectionSource{Float64}(50×25, RegularGridTree(51×26))"
            @test output == expected
        end
    end

    @testset "Output Conciseness Verification" begin

        @testset "ProjectionSource conciseness" begin
            x = range(-180, 180, length=721)
            y = range(-90, 90, length=361)
            data = rand(Float64, 720, 360)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            ps = SST.ProjectionSource(SST.RegularGridTree, ar)

            # Custom show method
            io_custom = IOBuffer()
            show(io_custom, MIME"text/plain"(), ps)
            custom_output = String(take!(io_custom))

            # Default show method
            io_default = IOBuffer()
            show(io_default, ps)
            default_output = String(take!(io_default))

            # Custom should be significantly more concise
            @test length(custom_output) < 100  # Should be under 100 characters
            @test length(custom_output) < length(default_output) / 5  # At least 5x more concise

            # Should contain essential information
            @test contains(custom_output, "Float64")
            @test contains(custom_output, "720×360")
            @test contains(custom_output, "RegularGridTree")
        end

        @testset "LazyProjectedDiskArray conciseness" begin
            x = range(0, 100, length=201)
            y = range(0, 50, length=101)
            data = rand(Float32, 200, 100)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            source = SST.ProjectionSource(SST.RegularGridTree, ar)

            target = SphericalSpatialTrees.ProjectionTarget(SST.ISEACircleTree, 5, 3)
            lpda = SphericalSpatialTrees.LazyProjectedDiskArray(source, target)

            # Custom show method
            io_custom = IOBuffer()
            show(io_custom, MIME"text/plain"(), lpda)
            custom_output = String(take!(io_custom))

            # Default show method
            io_default = IOBuffer()
            show(io_default, lpda)
            default_output = String(take!(io_default))

            # Custom should be significantly more concise
            @test length(custom_output) < 50  # Should be under 50 characters
            @test length(custom_output) < length(default_output) / 5  # At least 5x more concise

            # Should contain essential information
            @test contains(custom_output, "Float32")
            @test contains(custom_output, "32×32×10")
            @test contains(custom_output, "LazyProjectedDiskArray")
        end
    end

    @testset "Edge Cases and Error Handling" begin

        @testset "Very large arrays" begin
            # Test with very large dimensions (but don't actually create the data)
            x = range(-180, 180, length=10001)
            y = range(-90, 90, length=5001)
            # Create smaller actual data for memory efficiency
            data = rand(Float64, 100, 50)

            # Use smaller ranges for the actual array
            x_small = x[1:101]
            y_small = y[1:51]
            dims = (DD.X(x_small[1:end-1]), DD.Y(y_small[1:end-1]))
            ar = DD.DimArray(data, dims)
            ps = SST.ProjectionSource(SST.RegularGridTree, ar)

            # Should not throw exceptions
            io = IOBuffer()
            @test_nowarn show(io, MIME"text/plain"(), ps)
            output = String(take!(io))

            @test contains(output, "ProjectionSource{Float64}")
            @test contains(output, "100×50")
        end

        @testset "Different numeric types" begin
            # Test with various numeric types
            numeric_types = [Float16, Float32, Float64, Int8, Int16, Int32, Int64, UInt8, UInt16, UInt32]

            for T in numeric_types
                x = range(0, 1, length=11)
                y = range(0, 1, length=6)
                data = rand(T, 10, 5)

                dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
                ar = DD.DimArray(data, dims)
                ps = SST.ProjectionSource(SST.RegularGridTree, ar)

                io = IOBuffer()
                @test_nowarn show(io, MIME"text/plain"(), ps)
                output = String(take!(io))

                @test contains(output, "ProjectionSource{$T}")
                @test contains(output, "10×5")
            end
        end

        @testset "No exceptions during display" begin
            # Test that show methods never throw exceptions
            x = range(-180, 180, length=101)
            y = range(-90, 90, length=51)
            data = rand(Float64, 100, 50)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            source = SST.ProjectionSource(SST.RegularGridTree, ar)

            target = SST.ProjectionTarget(SST.ISEACircleTree, 4, 2)
            lpda = SST.LazyProjectedDiskArray(source, target)

            # Test all show methods don't throw
            @test_nowarn show(IOBuffer(), source.tree)
            @test_nowarn show(IOBuffer(), target.tree)
            @test_nowarn show(IOBuffer(), MIME"text/plain"(), source)
            @test_nowarn show(IOBuffer(), MIME"text/plain"(), lpda)
        end
    end

    @testset "Julia Array Display Convention Consistency" begin

        @testset "Dimension format consistency" begin
            # Test that dimensions are displayed with × separator like Julia arrays
            test_cases = [
                (10, 5, "10×5"),
                (100, 200, "100×200"),
                (1, 1, "1×1"),
                (1000, 500, "1000×500")
            ]

            for (dim1, dim2, expected_dims) in test_cases
                x = range(0, 1, length=dim1 + 1)
                y = range(0, 1, length=dim2 + 1)
                data = rand(Float64, dim1, dim2)

                dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
                ar = DD.DimArray(data, dims)
                ps = SST.ProjectionSource(SST.RegularGridTree, ar)

                io = IOBuffer()
                show(io, MIME"text/plain"(), ps)
                output = String(take!(io))

                @test contains(output, expected_dims)
            end
        end

        @testset "Type parameter format consistency" begin
            # Test that type parameters are in curly braces like Julia arrays
            x = range(0, 1, length=11)
            y = range(0, 1, length=6)
            data = rand(Float64, 10, 5)

            dims = (DD.X(x[1:end-1]), DD.Y(y[1:end-1]))
            ar = DD.DimArray(data, dims)
            ps = SST.ProjectionSource(SST.RegularGridTree, ar)

            io = IOBuffer()
            show(io, MIME"text/plain"(), ps)
            output = String(take!(io))
            # Should have type in curly braces
            @test contains(output, "{Float64}")
            @test startswith(output, "ProjectionSource{Float64}")
        end
    end
end