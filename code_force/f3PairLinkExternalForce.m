%RRR 2nd-Assur group dynamic statics

%input
%xe,ye-----position of E
%xsj,ysj---centroid sj of link j
%Frxe,Frye-----reaction force of kinematic pair upon E

%output
%Fcvtxj Fcvtyj ---- total reaction force converting E to centroid of j,cvt -- convert
%Mcvtfj ------ total reaction moment converting E to centroid of j

function [Fcvtxj,Fcvtyj,Mcvtfj] = f3PairLinkExternalForce(xe,ye,xsj,ysj,Frxe,Frye)

    Fcvtxj = - Frxe;
    Fcvtyj = - Frye;
    Mcvtfj = Frxe*(ye - ysj) - Frye*(xe - xsj);
    
end