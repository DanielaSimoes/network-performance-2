
N = 10; %número de simulações
results= zeros(1,N); %vetor com os N resultados de simulação

for it= 1:N
    results(it) = sim;
end

alfa= 0.1; %intervalo de confiança a 90%
media = mean(results);
termo = norminv(1-alfa/2)*sqrt(var(results)/N);
fprintf('resultado = %.2e +- %.2e\n',vpa(media),vpa(termo))

function [FinalResult] = sim()

    %Parameters initialization:
    N= 50;       % Number of mobile nodes
    W= 40;       % Radio range (in meters)
    S= 3;        % Maximum speed (in Km/h)
    delta= 1;    % Difference between consecutive time instants (in seconds)
    T= 3600;     % No. of time instants of the simulation
    %AP = [75 100;225 100];   % Coordinates of each AP
    %AP = [150 100];
    %AP = [50 100;150 100; 250 100];
    %AP = [75 50; 75 150; 225 50; 225 150];
    %AP = [50 75; 150 150; 250 75];

    %nAP = size(AP,1);    %Number of APs
    nAP = 1;
    S= S/3.6;            % Conversion of maximum speed to m/s
    results= zeros(1,T); % Initialization of the results array
    plotar = 0;  % if plotar = 1, node movement is visualized
                 % if plotar = 2, node movement and connectivity are visualized

    % Generation of initial coordinates of each mobile node position and speed:
    [pos,vel]= InitialRandom(N,S);
    h= waitbar(0,'Running simulation...');
    % Simulation cycle:
    for iter= 1:T
        waitbar(iter/T,h);
        % Update coordinates of each mobile node position and speed:
        [pos,vel]= UpdateCoordinates(pos,vel,delta);
        % Compute the node pairs with direct wireless links:
        L= ConnectedList(N,pos,W,AP);
        % Compute the percentage number of connected node pairs:
        results(iter)= AverageConnectedNodePairs(N,L,nAP);
        % Visualization of the simulation:
        if plotar>0
            visualize(pos,AP,L,plotar)
        end
    end
    delete(h)
    % Plot the simulation results:
    %figure(1)
    %plot(1:T,results)
    %axis([0 T 0 110])
    %xlabel('Time (seconds)');
    %ylabel('No. of connected nodes (%)')
    % Compute the final result: 
    FinalResult= mean(results);

end


function [pos,vel]= InitialRandom(N,S)
    %Computes a matrix �pos� of N rows and 2 columns with the coordinates of nodes (see Section 3.1) and a matrix �vel� of N rows and 2 columns with the initial horizontal and vertical components of the speed vector of each mobile node (see Section 3.2).
    pos= [50*randi([0 6],N/2,1) 200*rand(N/2,1)]; %primeira metade
    pos= [pos; 300*rand(N/2,1) 50*randi([0 4],N/2,1)]; %segunda metade
    vel_abs= S*rand(N,1);
    vel_angle= pi*randi([0 1],N/2,1) - pi/2; %primeira metade (de -pi/2 a pi/2)
    vel_angle= [vel_angle; pi*randi([0 1],N/2,1)]; %segunda metade (de 0 a pi)
    vel= [vel_abs.*cos(vel_angle) vel_abs.*sin(vel_angle)]; 
end

function [pos,vel]= UpdateCoordinates(pos,vel,delta)
    %Updates the matrices �pos� and �vel� based on their input values and delta.

    %atualizar os valores se a posição for menor que 0
    pos = pos + delta * vel;
    pos_x=pos(:,1);  %vamos buscar os valores da primeira coluna
    pos_y=pos(:,2);  %vamos buscar os valores da segunda coluna

    %mudar o sinal das velocidades cuja as posições saíam fora dos limites
    vel_x=vel(:,1);
    vel_y=vel(:,2);
    vel_x(pos_x<0 | pos_x>300)=-vel_x(pos_x<0 | pos_x>300);
    vel_y(pos_y<0 | pos_y>200)=-vel_y(pos_y<0 | pos_y>200);

    %devemos manter os valores dentro dos limites indicados
    pos_x(pos_x<0)=0; 
    pos_x(pos_x>300)=300;
    pos_y(pos_y<0)=0;
    pos_y(pos_y>200)=200;

    %atualizar a matriz com os valores corretos
    pos(:,1)=pos_x;
    vel(:,1)=vel_x;
    pos(:,2)=pos_y;
    vel(:,2)=vel_y;

end

function L= ConnectedList(N,pos,W,AP)
    %Computes a matrix �L� of 2 columns with the node pairs (mobile and AP nodes) such that their distance is not higher than W.
    % N nodes, Z APs
    
    W = W^2;
    L = [];
    nAPs=size(AP,1);
    pos = [pos; AP];
    
    for i=1:N
        for j=i+1:N+nAPs
            if ((pos(j,1)-pos(i,1))^2 + (pos(j,2)-pos(i,2))^2) < W
                L=[L; i j];
            end
        end
    end    
end

function o= AverageConnectedNodePairs(N,L,nAP)
%Computes a value �o� with the percentage number of connected node pairs based on the input matrix �L� of node pairs with direct links (see Section 4).

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
    
end

function visualize(pos,AP,L,plotar)
    N= size(pos,1);
    nAP= size(AP,1);
    plot(pos(1:N,1),pos(1:N,2),'o','MarkerEdgeColor','b','MarkerFaceColor','b')
    hold on
    plot(AP(1:nAP,1),AP(1:nAP,2),'s','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10)
    if plotar==2
        pos=[pos;AP];
        for i=1:size(L,1)
            plot([pos(L(i,1),1) pos(L(i,2),1)],[pos(L(i,1),2) pos(L(i,2),2)])
        end
    end
    axis([0 300 0 200])
    grid on
    set(gca,'xtick',0:50:300)
    set(gca,'ytick',0:50:200)
    drawnow
    hold off
end