//Kod zrodlowy programu symulujacego metapopulacje zlozona z subpopulacji zyjacych na rozlacznych platach powierzchni 

Uses Crt;
const
 ngen=14; //liczba genow
 seasons=5;
 k=13; //liczba powierzchni
 k0=4; //liczba powierzchni zrodlowych
 lambda=5;//6 mlodych
 ar=-0.1;
 br=0.61371;
 as=0.1;
 bs=-1.59453;
 pe=0.01;
 cc=0.00174;
 czas=10;//1000?
 lpowt=1;//1000? (czas=100,lpowt=100) -> symulacja trwala 6min na moim laptopie
 pmut=0.001; 
 skos=0.5; //pdb wydluzenia motywu. Dla skos = 0.5 rownie prawdopodobne wydluzenie co skrocenie.
 maxNAllel=1000;
 Ne=20;

Type stanosobnika=record
          NRsubpop:longint;
          NR:longint;
          plec: integer; //0=samica, 1=samiec
          geny: array[1..ngen,0..1] of longint; //pierwsza cyfra -> nr. genu; 0..1 -> 0 od matki, 1 od ojca; wartosc-> liczba powtorzen motywu 
          end;
 stanpopulacji=array[1..100000] of stanosobnika;
 nrSamcow=array[1..100000] of longint;
 tab = array[1..1000] of longint;
 liczebnosc=array[1..k] of longint;
 tablica=array[1..k, 1..k] of real;
 tablicaReal=array[1..k] of real;
 tablicaRealTime=array[1..czas,1..k] of real;
 GenDict=object
          labels: array[1..maxNAllel] of longint; //liczba powtorzen motywu w danym allelu
          values: array[1..maxNAllel] of longint; //liczba wystapien allelu
          max: longint;
          constructor init; 
          procedure add(lab:longint);
          procedure clear;
          procedure paste(dict:GenDict);
          function hasLabel(lab:longint):boolean;
          end;
 GenDict0=object
          labels: array[1..maxNAllel] of longint; //liczba powtorzen motywu w danym allelu
          values: array[1..maxNAllel] of real; //czestosc wystapien allelu
          max: longint;
          constructor init; 
          procedure add(lab:longint;val:real);
          procedure check(m,n:longint); 
          function getLabel(prob:real):longint;
          end;

 GenDictPop=array[1..k,1..ngen] of GenDict;
 GenDictPop0=array[1..k0,1..ngen] of GenDict0;


 constructor GenDict.init;
          begin 
          max:=0;
          end;
 procedure GenDict.add(lab:longint);
          var stop: boolean;
              i:longint;
          begin 
           stop:=false;
           for i:=1 to max do
           begin
            if labels[i]=lab then
             begin
             values[i]:=values[i]+1;
             stop:=true;
             break;
             end;
           end;
           if stop=false then
           begin
            labels[max+1]:=lab;
            values[max+1]:=1;
            max:=max+1;
           end;
          end;
 procedure GenDict.clear;
          begin
           max:=0;
          end;
 procedure GenDict.paste(dict:GenDict);
          var stop: boolean;
              i,j:longint;
          begin
           for i:=1 to dict.max do
           begin
            stop:=false; 
            for j:=1 to max do
            begin
             if dict.labels[i]=labels[j] then
             begin 
              values[j]:=values[j] + dict.values[i];
              stop:=true;
              break;
             end;
            end;
             if stop=false then 
             begin
              labels[max+1]:=dict.labels[i];
              values[max+1]:=dict.values[i];   
              max:=max+1;           
             end;
           end;
          end;
 function GenDict.hasLabel(lab:longint):boolean;
          var i:longint;
          begin
           hasLabel:=false;
           for i:=1 to max do 
            begin
             if labels[i]=lab then
              begin
               hasLabel:=true;
               break;
              end;
            end;
          end;

 constructor GenDict0.init;
          begin 
          max:=0;
          end;
 procedure GenDict0.add(lab:longint;val:real);
          begin
          max:=max+1;
          labels[max]:=lab;
          values[max]:=val;
          end;
 procedure GenDict0.check(m,n:longint); //Sprawdza, czy rozklad alleli jest poprawnie zdefiniowany
          var i,j:longint;
              sum:real;
          begin
          sum:=0;
          for i:=1 to max do
           begin
           for j:=i+1 to max do
            begin
            if labels[i]=labels[j] then writeln('allel ',labels[j],' w locus ',n, ' wystepuje dwa razy');
            end;
           sum:=sum+values[i];
           end;
          if sum<>1 then 
           begin
           if ((sum - 1)>0.1) or ((sum - 1)<-0.1) then writeln('suma prawdopodobienstw w populacji ',m,' i w locus ', n, ' wynosi ', sum:5:7, ' normuje prawdpopodobienstwa') ;
           for i:=1 to max do values[i]:=values[i]/sum;
           end;
          end;
 function GenDict0.getLabel(prob:real):longint;
          var i:longint;
              sum:real;
              stop:boolean;
          begin
           sum:=0;
           stop:=false;
           for i:=1 to max do
            begin
            sum:=sum+values[i];
            if prob<sum then 
             begin
             stop:=true;
             getLabel:=labels[i];
             break;
             end;
            end;
           if stop=false then
            begin
            writeln('blad w losowaniu! max=',max,' sum=',sum);
            getLabel:=0;
            end;
          end;
          


Var i, j, l, os, pot, t, powt, ll, dc, lm, s, Nt, locus, allel,N0,liczba,pairNumber: longint;
    N, Nm, Licz, Liczm, ost, lmigr : liczebnosc;
    pr, ps, sum, Ht, PopPairHt, freq, Fst: real; 
    polepow : tablicaReal;
    minodl,tabc : tablica;
    osob,potom,ojciec : stanosobnika;
    POP0,POP1,POPmigr : array[1..k] of stanpopulacji;
    samce0, samce1: array[1..k] of nrSamcow; 
    zrodla: array[1..k] of longint;
    znak : char;
    INFO, INFO1, INFO2, INFOsample, initgen, initpop, sampleSize : text; // INFOavg
    wiersz : string[k];
    word: string;
    allelsPop : GenDictPop;
    totalAllels,PopPairAllels: array[1..ngen] of GenDict;
    allels0: GenDictPop0;
    avgAllelRichness,avgPrivateAllels,He,Ho,Fis: tablicaReal;
    PopPairFst: array[1..round(k*(k-1)/2)] of real;
    //avgN, avgAAR, avgAPA, avgHe, avgHo, avgFis, avgFst: tablicaRealTime;
    //avgPPFst: array[1..czas,1..round(k*(k-1)/2)] of real;
    winter: boolean;

function normal(mi,sigma:real):real; //losuje liczbe z rozkladu normalnego: do losowania liczby potomstwa przy duzych lambda
Var alfa,r1,r2:real;
begin
  alfa:=2*pi*random;
  r1:=random;
  r2:=random;
  if r2<r1 then r2:=1-r2;
  normal:=sigma*cos(alfa)*sqrt(-2*ln(sqr(r2)))+mi;
end;

function combination(m,n:longint):tab; //losuje kombinacje bez powtorzen m liczb z posrod n (wynik w postaci array). Dla m>=n wybiera wszystkie (uwaga: liczba wynikow!).
 var i,j,l,r,count:longint;
     present:boolean;
 begin
 if m>=n then for i:=1 to n do combination[i]:=i else
  begin
  for i:=1 to m do
   begin
   r:=random(n-i+1)+1;
   writeln('random=',r);
   count:=0;
   j:=1;
   while(j<=r) do
    begin
    count:=count+1;
    present:=false;
    for l:=1 to m-1 do
     begin
     if combination[l]=count then present:=true;
     end;
    if present then else j:=j+1;
    end;
   combination[i]:=count; 
   end;
  end;
 end;
 
Function Lpotom(lambda:real):longint; //losuje liczbe potomkow z przesunietego o 1 rozkladu poissona o sredniej lambda
 var L:longint;
     x,y,LOT:real;
 begin
  if lambda<50 then
    begin
    L:=0;
    x:=exp(-lambda);
    y:=exp(-lambda);
    LOT:=random;
    while LOT>x do
      begin
      L:=L+1;
      y:=y*lambda/L;
      x:=x+y;
      end;
    end
  else L:=round(normal(lambda,sqrt(lambda)));
  LPOTOM:=L+1;
 end;


function POPDOC(j:longint; tabc:tablica):integer; //losuje populacje docelowa nie zalezy od k?
  var i:longint;
    czsum, los : real;
  begin
  los:=random;
  czsum:=0;
  i:=1;
  while los>czsum do
    begin
    czsum:=czsum+tabc[i,j];
    i:=i+1;
    end;
  POPDOC:=i-1;
  end;

function avgUnique(dictPop:GenDictPop): tablicaReal; //count averange number of private allels
 var i,j,l,n,lab:longint;
  stop:boolean;
 begin
 for i:=1 to k do //i -> nr populacji
  begin
   avgUnique[i]:=0;
   for j:=1 to ngen do //j -> nr genu
    begin 
     for l:=1 to dictPop[i][j].max do  //l -> nr allelu
      begin
       lab:=dictPop[i][j].labels[l];
       stop:=false;
       n:=0;
       while ((n < k) and (stop=false)) do //n -> nr populacji
        begin 
        n:=n+1;
        if ((dictPop[n][j].hasLabel(lab)) and (n<>i)) then
         begin
         stop:=true;
         break;
         end;
        end;
       if stop=false then avgUnique[i]:=avgUnique[i]+1;
      end;
    end;
   avgUnique[i]:=avgUnique[i]/ngen;
  end;
 end;

function expectHeterozigosity(dict:GenDict;N:longint): real;
 var l,sum:longint;
  homozigotes:real;
 begin
  sum:=0;
  if N=0 then expectHeterozigosity:=0 else
   begin
   homozigotes:=0;
   for l:=1 to dict.max do
    begin
    homozigotes:=homozigotes+(dict.values[l]*dict.values[l])/(4*N*N);
    sum:=sum + dict.values[l];
    end;
   if sum<>2*N then writeln ('sum = ',sum,' 2N = ',2*N);
   expectHeterozigosity:=1-homozigotes;
   end;
 end;

function expectPopHeterozigosity(dictPop:GenDictPop;N:liczebnosc): tablicaReal;
 var i,j:longint;
 begin
  for i:=1 to k do  
  begin
   expectPopHeterozigosity[i]:=0;
   for j:=1 to ngen do expectPopHeterozigosity[i]:=expectPopHeterozigosity[i]+expectHeterozigosity(dictPop[i][j],N[i]);
   expectPopHeterozigosity[i]:=expectPopHeterozigosity[i]/ngen;
  end;
 end;

{=========================================================================================================================
                                            B E G I N
==========================================================================================================================}

begin

{wczytanie danych o rozkladzie genotypow w populacjach zrodlowych z pliku INITgen.TXT}
 for i:=1 to k0 do for j:=1 to ngen do allels0[i,j].max:=0;
 assign(initgen,'INITgen.TXT');
 reset(initgen);
 //writeln('Wczytane czestosci alleli w populacjach zrodlowych');
 readln(initgen,word);
 //writeln(word);
 while not eof(initgen) do
  begin
  read(initgen,locus,znak,allel); 
  //write(locus,znak,allel);
  for i:=1 to k0 do
   begin
   read(initgen,znak,freq);
   //write(znak,freq:5:7);
   allels0[i,locus].add(allel,freq);
   end;
  readln(initgen);
  //writeln;
 end;
 for i:=1 to k0 do for j:=1 to ngen do allels0[i,j].check(i,j);
 close(initgen);  

{wczytanie danych o populacjach zrodlowych z pliku INITpop.TXT}
 assign(initpop,'INITpop.TXT');
 reset(initpop);
 readln(initpop,word);
 //writeln('Wczytane zrodel osobnikow w populacjach');
 for i:=1 to k do
  begin
  readln(initpop,liczba,znak,zrodla[i]);
  //writeln(liczba,znak,zrodla[i]);
  end;
 //writeln('koniec wczytywania');
 close(initpop);

{wczytanie danych z pliku POLEPOW.TXT}
  assign(INFO1,'POLEPOW.TXT');
  reset(INFO1);
  readln(INFO1,wiersz);
  for i:=1 to k do readln(INFO1,znak,POLEPOW[i]);
  close(INFO1);
  //writeln('wczytane pola powierzchni platow srodowiska');
  //for i:=1 to k do writeln(i,' ',polepow[i]:7:2,' ');

{wczytanie minimalnych odleglosci miedzy platami srodowiska}
  assign(INFO2,'MINODL.TXT');
  reset(INFO2);
  //writeln('wczytane odleglosci miedzy platami srodowiska');
  readln(INFO2,wiersz);
  for i:=1 to k do
    begin
    read(INFO2,znak);
    for j:=1 to k do 
     begin
     read(INFO2,MINODL[i,j]);
     //write(' ',MINODL[i,j]:7:2,' ');
     end;
    //writeln();
    readln(INFO2);
    end;
  close(INFO2);

{utworzenie tablicy tabc[i,j] na podstawie cc i minimalnych odleglosci}
  for j:=1 to k do
    begin
    sum:=0;
    for i:=1 to k do if j<>i then sum:=sum + exp(-cc*MINODL[i,j]); 
    for i:=1 to k do
      begin
      if (i=j) or (cc*MINODL[i,j]>300) then tabc[i,j]:=0 else tabc[i,j]:=exp(-cc*MINODL[i,j])/sum;
      end;
    end;

{utworzenie naglowka pliku wyjsciowego (zawierajacego wyniki poszczegolnych powtorzen symulacji)}
assign(INFO,'infodyn.txt');
 rewrite(INFO);
 append(INFO);
 write(INFO,'Powtorzenie czas ');
 for i:=1 to k do write(INFO,'N',i,' ');
 for i:=1 to k do write(INFO,'density',i,' ');
 for i:=1 to k do write(INFO,'avgAllelRichness',i,' ');
 for i:=1 to k do write(INFO,'avgPrivateAllels',i,' ');
 for i:=1 to k do write(INFO,'He',i,' ');
 for i:=1 to k do write(INFO,'Ho',i,' ');
 for i:=1 to k do write(INFO,'Fis',i,' ');
 for i:=1 to k do write(INFO,'Fst',i,' ');
 pairNumber:=0;
 for i:=1 to k do
  begin
   for j:=i+1 to k do 
    begin
    pairNumber:=pairNumber+1;
    write(INFO,'pairFst',i,'*',j,' ');
    end;
  end;
 writeln(INFO);
 close(INFO);

{utworzenie naglowka pliku wyjsciowego (zawierajacego wyniki srednich dla powtorzen)}
{assign(INFOavg,'infodynAvg.txt');
 rewrite(INFOavg);
 append(INFOavg);
 write(INFOavg,'czas ');
 for i:=1 to k do write(INFOavg,'N',i,' ');
 for i:=1 to k do write(INFOavg,'density',i,' ');
 for i:=1 to k do write(INFOavg,'avgAllelRichness',i,' ');
 for i:=1 to k do write(INFOavg,'avgPrivateAllels',i,' ');
 for i:=1 to k do write(INFOavg,'He',i,' ');
 for i:=1 to k do write(INFOavg,'Ho',i,' ');
 for i:=1 to k do write(INFOavg,'Fis',i,' ');
 for i:=1 to k do write(INFOavg,'Fst',i,' ');
 for pairNumber:=1 to round(k*(k-1)/2) do write(INFOavg,'pairFst',pairNumber,' ');
 writeln(INFOavg);
 close(INFOavg);}

{Utworzenie tablicy allelsPop}
 for i:=1 to k do
  begin
  for j:=1 to ngen do
   begin
   allelsPop[i,j].max:=0;
   for l:=1 to maxNAllel do
    begin
    allelsPop[i,j].labels[l]:=0;
    allelsPop[i,j].values[l]:=0;    
    end;
   end;
  end;

{Utworzenie tablicy totalAllels i PopPairAllels}
 for j:=1 to ngen do
  totalAllels[j].max:=0;
  PopPairAllels[j].max:=0;
  begin
  for l:=1 to maxNAllel do
   begin
   totalAllels[j].labels[l]:=0;
   totalAllels[j].values[l]:=0;
   PopPairAllels[j].labels[l]:=0;
   PopPairAllels[j].values[l]:=0;   
   end;
  end;

{===================================================================================================================================
						S T A R T  S Y M U L A C J I
====================================================================================================================================}
  randomize;

{wyzerowanie wskaznikow srednich dla wszystkich powtorzen}
 { for t:=1 to czas do
   begin
    for i:=1 to k do 
     begin
     avgN[t,i]:=0;
     avgAAR[t,i]:=0;
     avgAPA[t,i]:=0;
     avgHe[t,i]:=0;
     avgHo[t,i]:=0;
     avgFis[t,i]:=0;
     avgFst[t,i]:=0;
     end;
    for pairNumber:=1 to round(k*(k-1)/2) do avgPPFst[t,pairNumber]:=0;
   end;}
  


  for powt:=1 to lpowt do
    begin

    {wyzerowanie wskaznikow}    
    for j:=1 to ngen do totalAllels[j].clear;
    //for j:=1 to ngen do PopPairAllels[j].clear;
    Nt:=0;
    Ht:=0;
    Fst:=1;
    t:=0;
    for j:=1 to round(k*(k-1)/2) do PopPairFst[j]:=0;
    for i:=1 to k do 
      begin
      lmigr[i]:=0;
      avgAllelRichness[i]:=0;
      for j:=1 to ngen do allelsPop[i,j].clear;
      Nm[i]:=0;
      Ho[i]:=0;
      Fis[i]:=0;

      {utworzenie populacji poczatkowej}
     if zrodla[i]=0 then N0:=0 else N0:=round(0.1*Ne*polepow[i]);

      for os:=1 to N0 do
        begin
        osob.NRsubpop:=i;
        osob.nr:=os;
        osob.plec:=random(2);
	for j:=1 to ngen do
         begin
         for l:=0 to 1 do
          begin
          osob.geny[j,l]:=allels0[zrodla[i],j].getLabel(random);
          allelsPop[i,j].add(osob.geny[j,l]);
          end;
         if osob.geny[j,0]<>osob.geny[j,1] then Ho[i]:=Ho[i]+1;
         end;
        if osob.plec=1 then 
         begin
         Nm[i]:=Nm[i]+1;
         samce0[i][Nm[i]]:=os; //Zapisuje w tablicy nr samca w subpopulacji
         end;
        POP0[i][os]:=osob;
        end;
       ost[i]:=N0;
       N[i]:=N0;
      end;

    for i:=1 to k do 
     begin
     avgAllelRichness[i]:=0;
     for j:=1 to ngen do
      begin
      avgAllelRichness[i]:=avgAllelRichness[i]+allelsPop[i,j].max;
      totalAllels[j].paste(allelsPop[i,j]);
      end;
     avgAllelRichness[i]:=avgAllelRichness[i]/ngen;
     if N[i]<>0 then Ho[i]:=Ho[i]/(N[i]*ngen) else Ho[i]:=0;
     Nt:=Nt+N[i];
     end;
    avgPrivateAllels:=avgUnique(allelsPop); 
    He:=expectPopHeterozigosity(allelsPop,N);
    for j:=1 to ngen do Ht:=Ht+expectHeterozigosity(totalAllels[j],Nt);
    Ht:=Ht/ngen;

    for i:=1 to k do
     begin
     if He[i]<>0 then Fis[i]:=(He[i]-Ho[i])/He[i] else Fis[i]:=0;
     if Ht<>0 then Fst:=Fst - He[i]*N[i]/(Ht*Nt) else Fst:=0;
     end;

{wyliczanie wartosci Fst dla par populacji}
    pairNumber:=0;
    for i:=1 to k do
     begin
     for j:=i+1 to k do
      begin
      pairNumber:=pairNumber+1;
      PopPairHt:=0;
      for l:=1 to ngen do 
       begin
       PopPairAllels[l].clear;
       PopPairAllels[l].paste(allelsPop[i,l]);
       PopPairAllels[l].paste(allelsPop[j,l]);
       PopPairHt:=PopPairHt + expectHeterozigosity(PopPairAllels[l],N[i]+N[j])
       end;
      PopPairHt:=PopPairHt/ngen;
      if (N[i]+N[j])=0 then PopPairFst[pairNumber]:=1
      else if PopPairHt<>0 then PopPairFst[pairNumber]:=(PopPairHt-(He[i]*N[i]+He[j]*N[j])/(N[i]+N[j]))/PopPairHt 
      else PopPairFst[pairNumber]:=0;
      end;
     end;
      
    append(INFO);
    write(INFO,powt,' ',t,' ');
    for i:=1 to k do write(INFO,N[i],' ');
    for i:=1 to k do write(INFO,N[i]/polepow[i]:7:5,' ');
    for i:=1 to k do write(INFO,avgAllelRichness[i]:7:5,' ');
    for i:=1 to k do write(INFO,avgPrivateAllels[i]:7:5,' ');
    for i:=1 to k do write(INFO,He[i]:7:5,' ');
    for i:=1 to k do write(INFO,Ho[i]:7:5,' ');
    for i:=1 to k do write(INFO,Fis[i]:7:5,' ');
    write(INFO,Fst:7:5,' ');  
    for pairNumber:=1 to round(k*(k-1)/2) do write(INFO,PopPairFst[pairNumber]:7:5,' '); 
    writeln(INFO);
    ll:=0; //parametr przerywajacy symulacje, jesli ktoras z populacji jest zbyt liczna
    for i:=1 to k do if N[i]>100000 then ll:=1;
    close(INFO);

    {kolejne kroki czasowe}
    while ((t<czas)and(ll=0)) do
      begin
      t:=t+1;
      if ((t mod seasons = 4) or (t mod seasons = 0)) then winter:=true 
      else winter:=false;

      {wyzerowanie licznikow}
     Nt:=0;
     Ht:=0;
     Fst:=1;
     for j:=1 to ngen do totalAllels[j].clear;
      for i:=1 to k do
       begin 
       lmigr[i]:=0;
       avgAllelRichness[i]:=0;
       Ho[i]:=0;
       for j:=1 to ngen do allelsPop[i,j].clear;
       end;

      {przegladanie populacji}
      for i:=1 to k do
        begin
        Licz[i]:=0;
        Liczm[i]:=0;
        pr:=1/(1+exp(-(ar*N[i]/polepow[i]+br)));
        ps:=1/(1+exp(-(as*N[i]/polepow[i]+bs)));
        os:=0;
        while os<N[i] do
          begin
          os:=os+1;
          osob:=POP0[i][os];
          if ((not winter) and (osob.plec=0) and (random<pr) and (Nm[i]>0)) then for pot:=1 to Lpotom(lambda) do
            begin
            ojciec:=POP0[i][samce0[i][(random(Nm[i])+1)]]; //wybiera osobno ojca dla kazdego potomka
            potom.NRsubpop:=i;
            potom.nr:=ost[i]+1;
            potom.plec:=random(2);
            for j:=1 to ngen do
             begin
              potom.geny[j,0]:=osob.geny[j,(random(2))];
              potom.geny[j,1]:=ojciec.geny[j,(random(2))];
              for l:=0 to 1 do if random<pmut then
               begin
               if random<skos then potom.geny[j,l]:=potom.geny[j,l]+1
               else potom.geny[j,l]:=potom.geny[j,l]-1;
               end;
              allelsPop[i,j].add(potom.geny[j,0]);
              allelsPop[i,j].add(potom.geny[j,1]);
              if potom.geny[j,0]<>potom.geny[j,1] then Ho[i]:=Ho[i]+1;
             end; 
            ost[i]:=ost[i]+1;
            Licz[i]:=Licz[i]+1;
            POP1[i][Licz[i]]:=potom;
            if potom.plec=1 then
             begin
             Liczm[i]:=Liczm[i]+1;
             samce1[i][Liczm[i]]:=Licz[i];   
             end; 
            end;
          if Licz[i]=100000 then os:=N[i]+1 else
            begin
            if random<ps then else
              begin
              if random<pe then
                begin
                dc:=POPDOC(i, tabc);
                POPmigr[dc][lmigr[dc]+1]:=osob;
                lmigr[dc]:=lmigr[dc]+1;
                end
              else
                begin
                POP1[i][Licz[i]+1]:=osob;
                Licz[i]:=Licz[i]+1;
                for j:=1 to ngen do
                 begin
                 allelsPop[i,j].add(osob.geny[j,0]);
                 allelsPop[i,j].add(osob.geny[j,1]);
                 if osob.geny[j,0]<>osob.geny[j,1] then Ho[i]:=Ho[i]+1;
                 end;
                if osob.plec=1 then
                 begin
                 Liczm[i]:=Liczm[i]+1;
                 samce1[i][Liczm[i]]:=Licz[i];   
                 end;
                end;
              end;
            end;
          if Licz[i]=100000 then ll:=1;
          end;
        if Licz[i]=100000 then ll:=1;
        N[i]:=Licz[i];
        Nm[i]:=Liczm[i];
        end;
      for i:=1 to k do for lm:=1 to lmigr[i] do
        begin
        N[i]:=N[i]+1;
        POP1[i][N[i]]:=POPmigr[i][lm]; 
        osob:=POPmigr[i][lm];
        for j:=1 to ngen do
         begin
         allelsPop[i,j].add(osob.geny[j,0]);
         allelsPop[i,j].add(osob.geny[j,1]);
         if osob.geny[j,0]<>osob.geny[j,1] then Ho[i]:=Ho[i]+1;
         end;
        if osob.plec=1 then
         begin
         Nm[i]:=Nm[i]+1;
         samce1[i][Nm[i]]:=N[i];   
         end;
        end;

{statystyki}
      for i:=1 to k do 
       begin
       avgAllelRichness[i]:=0;
       for j:=1 to ngen do 
        begin
        avgAllelRichness[i]:=avgAllelRichness[i]+allelsPop[i,j].max;
        totalAllels[j].paste(allelsPop[i,j]);
        end;
       avgAllelRichness[i]:=avgAllelRichness[i]/ngen;
       if N[i]<>0 then Ho[i]:=Ho[i]/(N[i]*ngen) else Ho[i]:=0;
       Nt:=Nt+N[i];
       end;
      avgPrivateAllels:=avgUnique(allelsPop); 
      He:=expectPopHeterozigosity(allelsPop,N);
      for j:=1 to ngen do Ht:=Ht+expectHeterozigosity(totalAllels[j],Nt);
      Ht:=Ht/ngen;
      if Ht>1 then writeln('Ht!!!');
      if Ht<0 then writeln('-Ht!!!');
      for i:=1 to k do
       begin
       if He[i]<>0 then Fis[i]:=(He[i]-Ho[i])/He[i] else Fis[i]:=0;
       if ((Fis[i]<-1) or (Fis[i]>1)) then writeln('Fis[',i,']=',Fis[i]:5:7);
       if Ht*Nt<>0 then Fst:=Fst - He[i]*N[i]/(Ht*Nt) else Fst:=0;
       end;
      if Fst<0 then writeln('Fst<0 ','t=',t,' Ht=',Ht:5:7);

{wyliczanie wartości Fst dla par populacji}
      pairNumber:=0;
      for i:=1 to k do
       begin
       for j:=i+1 to k do
        begin
        pairNumber:=pairNumber+1;
        PopPairHt:=0;
        for l:=1 to ngen do 
         begin
         PopPairAllels[l].clear;
         PopPairAllels[l].paste(allelsPop[i,l]);
         PopPairAllels[l].paste(allelsPop[j,l]);
         PopPairHt:=PopPairHt + expectHeterozigosity(PopPairAllels[l],N[i]+N[j])
         end;
        PopPairHt:=PopPairHt/ngen;
        if (N[i]+N[j])=0 then PopPairFst[pairNumber]:=1
        else if PopPairHt<>0 then PopPairFst[pairNumber]:=(PopPairHt-(He[i]*N[i]+He[j]*N[j])/(N[i]+N[j]))/PopPairHt 
        else PopPairFst[pairNumber]:=0;
        end;
       end;

{sumy dla powtorzen}
      {for i:=1 to k do 
       begin
       avgN[t,i]:=avgN[t,i] + N[i];
       avgAAR[t,i]:=avgAAR[t,i] + avgAllelRichness[i];
       avgAPA[t,i]:=avgAPA[t,i] + avgPrivateAllels[i];
       avgHe[t,i]:=avgHe[t,i] + He[i];
       avgHo[t,i]:=avgHo[t,i] + Ho[i];
       avgFis[t,i]:=avgFis[t,i] + Fis[i];
       avgFst[t,i]:=avgFst[t,i] + Fst[i];
       end;
      for pairNumber:=1 to round(k*(k-1)/2) do avgPPFst[t,pairNumber]:=avgPPFst[t,pairNumber]+PopPairFst[pairNumber];}
      
      append(INFO);           
      write(INFO,powt,' ',t,' ');
      for i:=1 to k do write(INFO,N[i],' ');
      for i:=1 to k do write(INFO,N[i]/polepow[i]:7:5,' ');
      for i:=1 to k do write(INFO,avgAllelRichness[i]:7:5,' ');
      for i:=1 to k do write(INFO,avgPrivateAllels[i]:7:5,' ');
      for i:=1 to k do write(INFO,He[i]:7:5,' ');
      for i:=1 to k do write(INFO,Ho[i]:7:5,' ');
      for i:=1 to k do write(INFO,Fis[i]:7:5,' ');
      write(INFO,Fst:7:5,' ');
      for pairNumber:=1 to round(k*(k-1)/2) do write(INFO,PopPairFst[pairNumber]:7:5,' '); 
      writeln(INFO);
      for i:=1 to k do for os:=1 to N[i] do POP0[i][os]:=POP1[i][os];
      for i:=1 to k do for s:=1 to Nm[i] do samce0[i][s]:=samce1[i][s];
      gotoxy(1,whereY);
      if (t mod czas/100)=0 then write(powt,' ',t);
      end;
    end;
  close(INFO);

{srednie dla powtorzen}
 {append(INFOavg);
 for t:=1 to czas do
  begin
  for i:=1 to k do 
    begin
    avgN[t,i]:=avgN[t,i]/lpowt;
    avgAAR[t,i]:=avgAAR[t,i]/lpowt;
    avgAPA[t,i]:=avgAPA[t,i]/lpowt;
    avgHe[t,i]:=avgHe[t,i]/lpowt;
    avgHo[t,i]:=avgHo[t,i]/lpowt;
    avgFis[t,i]:=avgFis[t,i]/lpowt;
    avgFst[t,i]:=avgFst[t,i]/lpowt;
    end;
  for pairNumber:=1 to round(k*(k-1)/2) do avgPPFst[t,pairNumber]:=avgFst[t,pairNumber]/lpowt;
  write(INFOavg,t,' ');
  for i:=1 to k do write(INFOavg,avgN[t,i]:7:5,' ');
  for i:=1 to k do write(INFOavg,avgN[t,i]/polepow[i]:7:5,' ');
  for i:=1 to k do write(INFOavg,avgAAR[t,i]:7:5,' ');
  for i:=1 to k do write(INFOavg,avgAPA[t,i]:7:5,' ');
  for i:=1 to k do write(INFOavg,avgHe[t,i]:7:5,' ');
  for i:=1 to k do write(INFOavg,avgHo[t,i]:7:5,' ');
  for i:=1 to k do write(INFOavg,avgFis[t,i]:7:5,' ');
  for i:=1 to k do write(INFOavg,avgFst[t,i]:7:5,' ');
  for pairNumber:=1 to round(k*(k-1)/2) do write(INFOavg,avgPPFst[t,pairNumber]:7:5,' ');
  writeln(INFOavg);
 end;
 close(INFOavg);}
 writeln('Koniec');
 //repeat until Keypressed;
 //znak:=readkey;
end.

