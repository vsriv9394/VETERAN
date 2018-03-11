function flow_eval(scheme_conv,scheme_visc,Cells,SurfX,SurfY,SurfZ,gamma,Pr,mu)

    #---------------------------------------------------------------------------
    # Get the size for Cells vector
    #
    #   Structure of Cells vector
    #   [ρ, ρU, ρV, ρW, ρE, Cell_Volume]
    #---------------------------------------------------------------------------
    sCells  =   collect(size(Cells))

    #---------------------------------------------------------------------------
    # Create variable for output
    #---------------------------------------------------------------------------
    Flow    =   zeros(sCells[1],sCells[2],sCells[3],6)

    #---------------------------------------------------------------------------
    # Create dummy variables for left and right states
    #---------------------------------------------------------------------------
    stateL  =   zeros(6)
    stateR  =   zeros(6)

    #---------------------------------------------------------------------------
    # Create dummy variables for left and right gradients
    #---------------------------------------------------------------------------
    gradL   =   zeros(6,3)
    gradR   =   zeros(6,3)

    #---------------------------------------------------------------------------
    # Create dummy variable for surface between the left and right cells
    #---------------------------------------------------------------------------
    surf    =   [0.0,0.0,0.0,0.0]

    #Grads   =   reconstruct(Cells,SurfX,SurfY,SurfZ,gamma)

    
    ############################################################################
    # Calculate flows across ξ-surface
    ############################################################################
    
    for i=1:size(SurfX,1)
        for j=1:size(SurfX,2)
            for k=1:size(SurfX,3)

                #---------------------------------------------------------------
                # Set the boundary cell to the right by default
                #---------------------------------------------------------------
                ext_side        =   1

                #---------------------------------------------------------------
                # Set the boundary cell to the left for leftmost faces
                #---------------------------------------------------------------
                if i==1
                    ext_side    =   0
                end

                #---------------------------------------------------------------
                # Calculate the Area magnitude for the surface
                #---------------------------------------------------------------
                Area    =   sqrt(SurfX[i,j,k,1]^2+SurfX[i,j,k,2]^2+SurfX[i,j,k,3]^2)
                
                for m=1:5
                    #-----------------------------------------------------------
                    # Set the states to the left and right cell values
                    #-----------------------------------------------------------
                    stateL[m]   =   Cells[i  ,j+1,k+1,m]
                    stateR[m]   =   Cells[i+1,j+1,k+1,m]

                    #-----------------------------------------------------------
                    # Set the gradients to the left and right cell values
                    #-----------------------------------------------------------
                    #=
                    for mm=1:3
                        if m!=5
                            gradL[m,mm]  =   Grads[i  ,j+1,k+1,m,mm]
                            gradR[m,mm]  =   Grads[i+1,j+1,k+1,m,mm]
                        end
                    end
                    =#
                end

                #---------------------------------------------------------------
                # Set the surface area vector and surface boundary type
                #
                #   Structure of the SurfX vector
                #   [Area_x, Area_y, Area_z, Surface_Type]
                #---------------------------------------------------------------
                for m=1:4
                    surf[m]     =   SurfX[i,j,k,m]
                end

                #---------------------------------------------------------------
                # Evaluate net flow through the surface
                #---------------------------------------------------------------
                Flow_local  =   Area * setbcflux(scheme_conv,scheme_visc,gradL,gradR,stateL,stateR,gamma,Pr,mu,surf,ext_side)

                #---------------------------------------------------------------
                # Add and subtract fluxes to the cells right and left
                #---------------------------------------------------------------
                for m=1:5
                    
                    if i!=1

                        Flow[i  ,j+1,k+1,m] =   Flow[i  ,j+1,k+1,m]-Flow_local[m]/Cells[i  ,j+1,k+1,6]
                    
                    end

                    if i!=size(SurfX,1)
                        
                        Flow[i+1,j+1,k+1,m] =   Flow[i+1,j+1,k+1,m]+Flow_local[m]/Cells[i+1,j+1,k+1,6]
                    
                    end

                end

            end
        end
    end

    
    ############################################################################
    # Calculate flows across η-surface
    ############################################################################
    
    for i=1:size(SurfY,1)
        for j=1:size(SurfY,2)
            for k=1:size(SurfY,3)

                #---------------------------------------------------------------
                # Set the boundary cell to the right by default
                #---------------------------------------------------------------
                ext_side        =   1
                
                #---------------------------------------------------------------
                # Set the boundary cell to the left for leftmost faces
                #---------------------------------------------------------------
                if j==1
                    ext_side    =   0
                end

                #---------------------------------------------------------------
                # Calculate the Area magnitude for the surface
                #---------------------------------------------------------------
                Area    =   sqrt(SurfY[i,j,k,1]^2+SurfY[i,j,k,2]^2+SurfY[i,j,k,3]^2)
                
                for m=1:5
                    #-----------------------------------------------------------
                    # Set the states to the left and right cell values
                    #-----------------------------------------------------------
                    stateL[m]   =   Cells[i+1,j  ,k+1,m]
                    stateR[m]   =   Cells[i+1,j+1,k+1,m]

                    #-----------------------------------------------------------
                    # Set the gradients to the left and right cell values
                    #-----------------------------------------------------------
                    #=
                    for mm=1:3
                        if m!=5
                            gradL[m,mm]  =   Grads[i+1,j  ,k+1,m,mm]
                            gradR[m,mm]  =   Grads[i+1,j+1,k+1,m,mm]
                        end
                    end
                    =#
                end

                #---------------------------------------------------------------
                # Set the surface area vector and surface boundary type
                #
                #   Structure of the SurfX vector
                #   [Area_x, Area_y, Area_z, Surface_Type]
                #---------------------------------------------------------------
                for m=1:4
                    surf[m]     =   SurfY[i,j,k,m]
                end

                #---------------------------------------------------------------
                # Evaluate net flow through the surface
                #---------------------------------------------------------------
                Flow_local  =   Area * setbcflux(scheme_conv,scheme_visc,gradL,gradR,stateL,stateR,gamma,Pr,mu,surf,ext_side)

                #---------------------------------------------------------------
                # Add and subtract fluxes to the cells right and left
                #---------------------------------------------------------------
                for m=1:5

                    if j!=1

                        Flow[i+1,j  ,k+1,m] =   Flow[i+1,j  ,k+1,m]-Flow_local[m]/Cells[i+1,j  ,k+1,6]

                    end
                    if j!=size(SurfY,2)

                        Flow[i+1,j+1,k+1,m] =   Flow[i+1,j+1,k+1,m]+Flow_local[m]/Cells[i+1,j+1,k+1,6]

                    end

                end

            end
        end
    end

    
    ############################################################################
    # Calculate flows across ζ-surface
    ############################################################################
    
    for i=1:size(SurfZ,1)
        for j=1:size(SurfZ,2)
            for k=1:size(SurfZ,3)

                #---------------------------------------------------------------
                # Set the boundary cell to the right by default
                #---------------------------------------------------------------
                ext_side        =   1
                
                #---------------------------------------------------------------
                # Set the boundary cell to the left for leftmost faces
                #---------------------------------------------------------------
                if k==1
                    ext_side    =   0
                end

                #---------------------------------------------------------------
                # Calculate the Area magnitude for the surface
                #---------------------------------------------------------------
                Area    =   sqrt(SurfZ[i,j,k,1]^2+SurfZ[i,j,k,2]^2+SurfZ[i,j,k,3]^2)
                
                for m=1:5
                    #-----------------------------------------------------------
                    # Set the states to the left and right cell values
                    #-----------------------------------------------------------
                    stateL[m]   =   Cells[i+1,j+1,k  ,m]
                    stateR[m]   =   Cells[i+1,j+1,k+1,m]

                    #-----------------------------------------------------------
                    # Set the gradients to the left and right cell values
                    #-----------------------------------------------------------
                    #=
                    for mm=1:3
                        if m!=5
                            gradL[m,mm]  =   Grads[i+1,j+1,k  ,m,mm]
                            gradR[m,mm]  =   Grads[i+1,j+1,k+1,m,mm]
                        end
                    end
                    =#
                end

                #---------------------------------------------------------------
                # Set the surface area vector and surface boundary type
                #
                #   Structure of the SurfX vector
                #   [Area_x, Area_y, Area_z, Surface_Type]
                #---------------------------------------------------------------
                for m=1:4
                    surf[m]     =   SurfZ[i,j,k,m]
                end

                #---------------------------------------------------------------
                # Evaluate net flow through the surface
                #---------------------------------------------------------------
                Flow_local  =   Area * setbcflux(scheme_conv,scheme_visc,gradL,gradR,stateL,stateR,gamma,Pr,mu,surf,ext_side)

                #---------------------------------------------------------------
                # Add and subtract fluxes to the cells right and left
                #---------------------------------------------------------------
                for m=1:5

                    if k!=1

                        Flow[i+1,j+1,k  ,m] =   Flow[i+1,j+1,k  ,m]-Flow_local[m]/Cells[i+1,j+1,k  ,6]

                    end

                    if k!=size(SurfZ,3)

                        Flow[i+1,j+1,k+1,m] =   Flow[i+1,j+1,k+1,m]+Flow_local[m]/Cells[i+1,j+1,k+1,6]

                    end

                end

            end
        end
    end

    return Flow

end
