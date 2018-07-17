function o= AverageConnectedNodePairs(N,L,nAP)
    repeat = 1;
    
    labels_nodes = zeros(N+nAP,1); %tudo com zeros
    
    for i=1:nAP
        labels_nodes(N+i,1) = 1; %colocar os APs com uns
    end
    
    while repeat == 1
        repeat = 0;
        for i=1:size(L,1) %percorrer todos os pares de nós     
            if(labels_nodes(L(i,1)) ~= labels_nodes(L(i,2)))
                %verificar se as posições estão a 0, se estiverem, alterar
                %para 1, APS têm semplre o valor 1
                if labels_nodes(L(i,1)) == 0 
                    labels_nodes(L(i,1)) = 1;
                else
                    labels_nodes(L(i,2)) = 1;
                end
                %repeat fica igual a 1
                repeat = 1;
            end
        end
    end
    
    %verificar quantos valores estão a 1 
    count = sum(labels_nodes(1:N));
    %fazer o calculo do O
    o = (count/N)*100;
    