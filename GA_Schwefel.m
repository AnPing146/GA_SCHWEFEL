clc;clear;close all;

ITERATION = 20;
REPEAT = 10;
POPULATION = 96;

chromosome = zeros(POPULATION, 2, ITERATION);
result = zeros(REPEAT, ITERATION, 3);       % 3-type average

%% fitness function
x1 = -500:10:500;
x2 = -500:10:500;
y = zeros(length(x1), length(x2));

ii=1;
for i = x1
    jj=1;
    for j = x2
        y(ii, jj) = SCHWEFEL(i, j);
        jj=jj+1;
    end
    ii=ii+1;
end
%% file read (result)
fptr = fopen("C:\Users\71713\source\repos\Optima\GA_midterm\GA_midterm\GA_result.txt",'r');

for i=1:3
    fgetl(fptr);

    for j=1:REPEAT
        for k=1:ITERATION
            result(j,k,i) = fscanf(fptr,'%f ',[1,1]);
        end
    end
end

fclose(fptr);

%plot result
for i=1:3
    figure(i);

    for j=1:REPEAT
        plot(1:20,result(j, :, i));hold on;

        xlabel('iteration');ylabel('fitness');
        switch i
            case 1
                title("total average");
            case 2
                title("top-48 average");
            case 3
                title("top-20 average");
        end
    end
end


%% file read (chromosome)
fptr = fopen('C:\Users\71713\source\repos\Optima\GA_midterm\GA_midterm\GA_chromosome.txt','r');

for i=1:ITERATION
    fgetl(fptr);

    for j=1:POPULATION
        chromosome(j,:,i) = fscanf(fptr,'(%f,%f) ',[1, 2]);
    end
end

fclose(fptr);

% plot chromosome
y2 = zeros(POPULATION, 1, ITERATION);

for i=1:ITERATION
    for j=1:POPULATION
        y2(j, 1, i) = SCHWEFEL(chromosome(j, 1, i), chromosome(j, 2,i));
    end
end

f = figure(4);f.Position=[50 50 800 600];
pause(2);
for i = 1:ITERATION
    surf(x1, x2, y); hold on;                                          %plot fitness function
    title("iteration:"+i);colorbar;xlabel('x1');ylabel('x2');
    g = plot3(chromosome(:,1,i), chromosome(:,2,i), y2(:,1,i), 'r*');  %plot chromosome
   
    view(0, 90);

    pause(0.1);

    if(i<ITERATION)
        delete(g);
    end
end

