fprintf("\nUpdate these Variables: (if none, YAY!)\n")
if round(CD_0,4) ~= round(CD_0i,4)
    fprintf("Update CD_0i to %.4f\n", CD_0)
end

if round(CL_R,4) ~= round(CL_Ri,4)
    fprintf("Update CL_Ri to %.4f\n", CL_R)
end

if round(etaP_CR,4) ~= round(etaP_CRi,4)
    fprintf("Update etaP_CRi to %.4f\n", etaP_CR)
end

if round(etaP_TO,4) ~= round(etaP_TOi,4)
    fprintf("Update etaP_TOi to %.4f\n", etaP_TO)
end

if round(etaP_M,4) ~= round(etaP_Mi,4)
    fprintf("Update etaP_Mi to %.4f\n", etaP_M)
end

if round(etaP_CL,4) ~= round(etaP_CLi,4)
    fprintf("Update etaP_CLi to %.4f\n", etaP_CL)
end
