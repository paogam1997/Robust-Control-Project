%% Funzione per trovare l'upper bound di un insieme di diagrammi di Bode
function [W]=upperbound(fun_trasf,ordine,numero_di_punti,ord_rel)
    G=fun_trasf;
    ord=ordine;
    n=numero_di_punti;
    ord_rel=[];
    if  isempty(ord_rel)==1
        RD=0;
    else
        RD=ord_rel;
    end
    figure;
    bode(G);

    
    
    [freq,resp_db]=ginput(n); 
    %for i=1:n
        %resp(i)=10^(resp_db(i)/20); 
    %end
    resp=10.^(resp_db/20);
    sys=frd(resp,freq);  
    W=fitmagfrd(sys,ord,RD); 
    W=tf(W)
    hold on
    bode(W);
    hold off
end