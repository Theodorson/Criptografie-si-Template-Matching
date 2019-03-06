#include <stdio.h>
#include <stdlib.h>
#include <math.h>

///STRUCTURA FOLOSITA LA XORAREA PIXELILOR CU UN INTREG FARA SEMN
typedef union
{
    unsigned int x;
    char octet[4];
} octeti;

///STRUCTURA FOLOSITA LA MATRICEA DE TIP PIXEL
typedef struct
{
    unsigned char octeti[3];
} pixel;


///STRUCTURA FOLOSITA LA RETINEREA DETECTIILOR
typedef struct
{
    double cor;
    int coordonate[2];
    int cif;
    int ok;
} detectii;

///STRUCTURA FOLOSITA LA COLORAREA CHENARELOR
typedef struct
{
    unsigned char culoare[3];
} contur;

///FUNCTIA CMPDESCRESCATOR FOLOSITA IN QSORT
int cmpDescrescator (const void *a,const void *b)
{
    detectii va=*(detectii*)a;
    detectii vb=*(detectii*)b;
    if (va.cor<vb.cor) return 1;
    if (va.cor>vb.cor) return -1;
    return 0;
}

unsigned int* generator_XORSHIFT32(int x,int W,int H)
{
    unsigned int r,k,seed=x;
    unsigned int *aleator=(unsigned int*)malloc((2*W*H-1)*sizeof(unsigned int));
    r=seed;
    for (k=0; k<2*W*H-1; k++)
    {
        r=r^r<<13;
        r=r^r>>17;
        r=r^r<<5;
        aleator[k]=r;
    }
    return aleator;
}

unsigned char* incarcare_im(char *cale_imagine) ///LINIARIZARE IMAGINE
{
    FILE *in=fopen(cale_imagine,"rb");
    if (in==NULL) printf("EROARE");
    int w,h,*W,*H,padding=0,p=0,j;
    dimensiuni_imagine(cale_imagine,&W,&H);
    unsigned char *incarcare=NULL,*aux=NULL;
    int n=0,i=0;
    unsigned char x;
    fseek(in,54,SEEK_SET);
    w=W;
    h=H;
    if (w%4!=0)
    padding=4-(3*w)%4;
    else
    padding=0;

        i=0;
        p=0;
        for (j=0; j<h; j++)
        {
            while (fread(&x,sizeof(unsigned char),1,in)==1)
            {
                if (i!=3*w)
                {
                    n++;
                    i++;
                    aux=(unsigned char*)realloc(incarcare,n*sizeof(unsigned char));
                    if (aux==NULL)
                    {
                        free(incarcare);
                        return 0;
                    }
                    else
                    {
                        incarcare=aux;
                        incarcare[n-1]=x;
                    }
                }
                else
                {   if (padding!=0)
                   {
                     while (p<padding-1)
                    {
                        fread(&x,sizeof(unsigned char),1,in);
                        p++;
                    }
                    p=0;
                    i=0;
                   }
                   else
                   {
                       n++;
                       i++;
                       aux=(unsigned char*)realloc(incarcare,n*sizeof(unsigned char));
                    if (aux==NULL)
                    {
                        free(incarcare);
                        return 0;
                    }
                    else
                    {
                        incarcare=aux;
                        incarcare[n-1]=x;
                    }
                   }
                }
            }
        }
    fclose(in);
    return incarcare;
}


void creare_fisier(char *cale_imagine,char *cale_imagine_criptata, int W, int H,unsigned char *pixeli) ///CREARE IMAGINE CRIPTATA/DECRIPTATA
{
    unsigned char *header=Header_im(cale_imagine),x;
    FILE *out=fopen(cale_imagine_criptata,"wb"),*in=fopen(cale_imagine,"rb");
    if (out==NULL)
    {
        printf("EROARE");
        return 0;
    }
    if(in==NULL)
    {
        printf("EROARE");
        return 0;
    }
    int i=0,j,nrpadding;
    while (i<54)
    {
        fwrite(&header[i],sizeof(unsigned char),1,out);
        i++;
    }

    if (W%4!=0)
        nrpadding=4-(3*W)%4;
    else
        nrpadding=0;
    i=0,j=0;
    while (i<=3*W*H)
    {

        if (j==3*W && nrpadding!=0)
        {
            j=0;
            while (j<nrpadding)
            {
                fread(&x,sizeof(unsigned char),1,in);
                fwrite(&x,sizeof(unsigned char),1,out);
                j++;
            }
            j=0;
        }
        if (i!=3*W*H)
            fwrite(&pixeli[i],sizeof(unsigned char),1,out);

        fread(&x,sizeof(unsigned char),1,in);
        j++;
        i++;
    }
    fclose(out);
    fclose(in);
    free(header);
}

unsigned int* alg_Durstenfeld(char *cale_cheie,int W,int H) /// FUNCTIE CARE RETURNEAZA O PERMUTARE ALEATOARE
{
    unsigned int x;
    FILE *in=fopen(cale_cheie,"r");
    fscanf(in,"%d",&x);
    fclose(in);
    unsigned int *aleator=generator_XORSHIFT32(x,W,H);
    unsigned int r,aux;
    unsigned int *permutare=(unsigned int *)malloc(W*H*sizeof(unsigned int));
    int k;
    for (k=0; k<W*H; k++)
    {
        permutare[k]=k;
    }
    for (k=W*H-1; k>=1; k--)
    {
        r=aleator[k]%(k+1);
        aux=permutare[r];
        permutare[r]=permutare[k];
        permutare[k]=aux;
    }
    return permutare;
}

unsigned char* alg_criptare(char *cale_imagine,char *cale_cheie,char *cale_imagaine_criptata,int W,int H)
{
   ///FUNCTIE CARE CRIPTEAZA SI RETURNEAZA UN VECTOR CU PIXELII CRIPTATI
    unsigned int x;
    octeti SV;
    FILE *in=fopen(cale_cheie,"r");
    fscanf(in,"%lu",&x);
    fscanf(in,"%lu",&SV.x);
    fclose(in);
    unsigned int *aleator=generator_XORSHIFT32(x,W,H),*permutare_initiala=alg_Durstenfeld(cale_cheie,W,H);
    unsigned char *pixeli_perm=incarcare_im(cale_imagine),*pixeli_criptati=incarcare_im(cale_imagine);
    int pozitie_pixel=0,j,i,k,aux;


    /// PERMUTAREA PIXELILOR
    while (pozitie_pixel<=W*H-1)
    {

        pixeli_perm[3*permutare_initiala[pozitie_pixel]+0]=pixeli_criptati[3*pozitie_pixel+0];
        pixeli_perm[3*permutare_initiala[pozitie_pixel]+1]=pixeli_criptati[3*pozitie_pixel+1];
        pixeli_perm[3*permutare_initiala[pozitie_pixel]+2]=pixeli_criptati[3*pozitie_pixel+2];
        pozitie_pixel++;
    }
/// SCHIMBAREA VALORILOR PIXELILOR
    pozitie_pixel=0;
    octeti aux_xor;
    while (pozitie_pixel<=W*H-1)
    {
        for (j=0; j<3; j++)
        {
            if (pozitie_pixel==0)
            {
                aux_xor.x=aleator[W*H];
                pixeli_criptati[3*pozitie_pixel+j]=SV.octet[j]^pixeli_perm[3*pozitie_pixel+j]^aux_xor.octet[j];
            }
            else
            {
                aux_xor.x=aleator[W*H+pozitie_pixel];
                pixeli_criptati[3*pozitie_pixel+j]=pixeli_perm[3*pozitie_pixel+j]^pixeli_criptati[3*(pozitie_pixel-1)+j]^aux_xor.octet[j];
            }
        }
        pozitie_pixel++;
    }
    free(aleator);
    free(permutare_initiala);
    return pixeli_criptati;
}

unsigned char* alg_decriptare (char *cale_imagine,char *cale_cheie,char *cale_imagine_criptata,int W,int H,unsigned char *pixeli_criptati)
{
    ///FUNCTIE CARE DECRIPTEAZA SI RETURNEAZA UN VECTOR CU PIXELII DECRIPTATI
    octeti SV;
    unsigned int x;
    FILE *cheie=fopen(cale_cheie,"r");
    if (cheie==NULL)
    {
        printf("EROARE");
        return 0;
    }
    fscanf(cheie,"%lu",&x);
    fscanf(cheie,"%lu",&SV.x);
    fclose(cheie);
    unsigned int *permutare_initiala=alg_Durstenfeld(cale_cheie,W,H),*aleator=generator_XORSHIFT32(x,W,H),*permutare_inversa=alg_Durstenfeld(cale_cheie,W,H);
    unsigned char *pixeli_decriptati1=alg_criptare(cale_imagine,cale_cheie,cale_imagine_criptata,W,H),*pixeli_decriptati2=alg_criptare(cale_imagine,cale_cheie,cale_imagine_criptata,W,H);
    int j=0,i=0;

    ///INVERSA
    for (i=0; i<W*H; i++)
        permutare_inversa[permutare_initiala[i]]=i;

    octeti aux_xor;
    int pozitie_pixel=0;

    //SCHIMBAREA VALORILOR PIXELILOR
    while (pozitie_pixel<=W*H-1)
    {
        for (j=0; j<3; j++)
        {
            if (pozitie_pixel==0)
            {
                aux_xor.x=aleator[W*H];
                pixeli_decriptati1[3*pozitie_pixel+j]=SV.octet[j]^pixeli_criptati[3*pozitie_pixel+j]^aux_xor.octet[j];
            }
            else
            {
                aux_xor.x=aleator[W*H+pozitie_pixel];
                pixeli_decriptati1[3*pozitie_pixel+j]=pixeli_criptati[3*pozitie_pixel+j]^pixeli_criptati[3*(pozitie_pixel-1)+j]^aux_xor.octet[j];
            }
        }
        pozitie_pixel++;
    }
    pozitie_pixel=0;

    //PERMUTAREA PIXELILOR
    while (pozitie_pixel<=W*H-1)
    {
        pixeli_decriptati2[3*permutare_inversa[pozitie_pixel]+0]=pixeli_decriptati1[3*pozitie_pixel+0];
        pixeli_decriptati2[3*permutare_inversa[pozitie_pixel]+1]=pixeli_decriptati1[3*pozitie_pixel+1];
        pixeli_decriptati2[3*permutare_inversa[pozitie_pixel]+2]=pixeli_decriptati1[3*pozitie_pixel+2];
        pozitie_pixel++;
    }
    free(pixeli_decriptati1);
    free(permutare_initiala);
    free(permutare_inversa);
    return pixeli_decriptati2;
}

void hi_2(char *cale_imagine,int W,int H)
{   ///TESTUL CHI-PATRAT
    int i,j;
    float s=0;
    unsigned char *continut=incarcare_im(cale_imagine);
    int k;
    float a=256;
    float f=(W*H)/a;
    for(i=0; i<3; i++)
    {
        s=0;
        float frecventa[256]={0};
        for (j=i; j<3*W*H; j=j+3)
        {
            frecventa[continut[j]]++;
        }

        for (k=0; k<=255; k++)
        {
            s=s+(pow((frecventa[k]-f),2)/f);
        }
        if (i==0)
        printf("\nR=%.2f \n",s);
        if (i==1)
        printf("G=%.2f \n",s);
        if (i==2)
        printf("B=%.2f \n",s);
    }
    printf("\n");
    free(continut);
}

void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie)
{
    FILE *fin, *fout;
    unsigned int dim_img, latime_img, inaltime_img;
    unsigned char pRGB[3], header[54], aux;

    //printf("nume_fisier_sursa = %s \n",nume_fisier_sursa);

    fin = fopen(nume_fisier_sursa, "rb");
    if(fin == NULL)
    {
        //printf("nu am gasit imaginea sursa din care citesc");
        return;
    }

    fout = fopen(nume_fisier_destinatie, "wb+");

    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);
    //printf("Dimensiunea imaginii in octeti: %u\n", dim_img);

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);
   // printf("Dimensiunea imaginii in pixeli (latime x inaltime): %u x %u\n",latime_img, inaltime_img);

    //copiaza octet cu octet imaginea initiala in cea noua
    fseek(fin,0,SEEK_SET);
    unsigned char c;
    while(fread(&c,1,1,fin)==1)
    {
        fwrite(&c,1,1,fout);
        fflush(fout);
    }
    fclose(fin);

    //calculam padding-ul pentru o linie
    int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    //printf("padding = %d \n",padding);

    fseek(fout, 54, SEEK_SET);
    int i,j;
    for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            //citesc culorile pixelului
            fread(pRGB, 3, 1, fout);
            //fac conversia in pixel gri
            aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
            pRGB[0] = pRGB[1] = pRGB[2] = aux;
            fseek(fout, -3, SEEK_CUR);
            fwrite(pRGB, 3, 1, fout);
            fflush(fout);
        }
        fseek(fout,padding,SEEK_CUR);
    }
    fclose(fout);
}

void dimensiuni_imagine(char *cale_imagine, int *W, int *H)
{

    FILE *in=fopen(cale_imagine,"rb");
    if (in==NULL) printf("EROARE");
    int i=0;
    unsigned char a;
    while (i<=54 && fread(&a,sizeof(unsigned char),1,in)==1)
    {
        if (i==17)
        {
            fread(&(*W),sizeof(int),1,in);
            i=i+4;
        }
        if (i==21)
        {
            fread(&(*H),sizeof(int),1,in);
            i=i+4;
        }
        if (i!=17 && i!=21) i++;
    }
    fclose(in);
}

int Header_im(char *nume_imagine)
{
    FILE *in=fopen(nume_imagine,"rb");
    if (in==NULL) printf("EROARE");
    unsigned char *header=NULL;
    int i=0;
    header=(unsigned char*)malloc(54*sizeof(unsigned char));
    while(i<54)
    {
        fread(&header[i],sizeof(unsigned char),1,in);
        i++;
    }
    fclose(in);
    return header;
}


pixel** img_matrice(char* nume_img,int W, int H)
{
    FILE *in=fopen(nume_img,"rb");
    if (in==NULL)
    {
        printf("EROARE");
        return 0;
    }
    int i,j,padding,k;
    pixel **matrice_img=(pixel**)malloc(H*sizeof(pixel*));
    unsigned char x;
    for (i=0; i<H; i++)
        matrice_img[i]=(pixel*)malloc(W*sizeof(pixel));
    if (W%4!=0)
        padding=4-(3*W)%4;
    else
        padding=0;
    fseek(in,54,SEEK_SET);
    ///TRANSFORMAREA IMAGINEI INTR-O MATRICE DE TIP PIXEL
    for (i=0; i<H; i++)
    {
        for (j=0; j<W; j++)
        {
            fread(&matrice_img[i][j],sizeof(pixel),1,in);
        }
        for (j=0; j<padding; j++)
            fread(&x,sizeof(unsigned char),1,in);
    }
    fclose (in);
    return matrice_img;
}

void copie_imagine(char *nume_imagine,char *nume_imagine_modificata,int W, int H)
{
    pixel **continut=img_matrice(nume_imagine,W,H);
    unsigned char *header=Header_im(nume_imagine);
    int i=0,j=0,nrpadding,a=0;
    FILE *out=fopen(nume_imagine_modificata,"wb");
    if (out==NULL)
    {
        printf("EROARE");
        return 0;
    }
    while (i<54)
    {
        fwrite(&header[i],sizeof(unsigned char),1,out);
        i++;
    }
    i=0;
    if (W%4!=0)
        nrpadding=4-(3*W)%4;
    else
        nrpadding=0;
    i=0,j=0;
    for (i=0; i<H; i++)
    {
        for (j=0; j<W; j++)
            fwrite(&continut[i][j],sizeof(pixel),1,out);
        for (j=0; j<nrpadding; j++)
            fwrite(&a,sizeof(int),1,out);
    }
    fclose(out);
    for (i=0; i<H; i++)
        free(continut[i]);
    free(continut);
}

void template_matching(char* nume_img,char* nume_sablon,int W,int H,int w,int h,float ps,int NRX,contur C[],detectii **D,int *nr)
{
    pixel **fereastra=img_matrice(nume_img,W,H);
    pixel **sablon=img_matrice(nume_sablon,w,h);
    unsigned int i,j,x,y,k,a,p,numarare;
    double S_medie=0,corelatie=0,f_medie=0,n=w*h,Suma=0,deviatie_sablon=0,deviatie_fereastra=0,S=0,f_i=0;
    detectii *aux=NULL;
    ///CALCULAREA MEDIEI PIXELILOR SABLONULUI
    for (i=0; i<h; i++)
        for (j=0; j<w; j++)
        {
            S=sablon[i][j].octeti[1];
            S_medie=S_medie+S;
        }
    S_medie=S_medie/n;

///CALCULAREA DEVIATIEI SABLONULUI
    for (i=0; i<h; i++)
        for (j=0; j<w; j++)
        {
            S=sablon[i][j].octeti[1];
            deviatie_sablon=deviatie_sablon+pow((S-S_medie),2);
        }
    deviatie_sablon=sqrt((1/(n-1))*deviatie_sablon);



    x=0;
    y=0;
    numarare=0;
    while (x<H)
    {
        while (y<W)
        {
            ///CALCULAREA MEDIEI PIXELILOR IN FEREASTRA I
            for (i=x; i<h+x; i++)
                for (j=y; j<w+y; j++)
                {
                    if (i>=H || j>=W)
                    {
                        f_i=0;
                    }
                    else
                    {
                        f_i=fereastra[i][j].octeti[1];
                    }
                    f_medie=f_medie+f_i;
                }
            f_medie=f_medie/n;
            ///CALCULAEA DEVIATIEI FERESTREI
            for (i=x; i<h+x; i++)
                for (j=y; j<w+y; j++)
                {
                    if (i>=H || j>=W)
                    {
                        f_i=0;
                    }
                    else
                    {
                        f_i=fereastra[i][j].octeti[1];
                    }
                    deviatie_fereastra=deviatie_fereastra+pow((f_i-f_medie),2);
                }
            deviatie_fereastra=sqrt((1/(n-1))*deviatie_fereastra);
            i=0;
            ///APLIC FORMULA DIN ENUNT
            for (k=x; k<h+x; k++)
            {
                j=0;
                for (p=y; p<w+y; p++)
                {
                    if (k>=H || p>=W)
                    {
                        f_i=0;
                    }
                    else
                    {
                        f_i=fereastra[k][p].octeti[1];
                    }
                    S=sablon[i][j].octeti[1];
                    Suma=Suma+((1/(deviatie_sablon*deviatie_fereastra))*(f_i-f_medie)*(S-S_medie));
                    j++;
                }
                i++;
            }
            corelatie=(1/n)*Suma*4;
            Suma=0;
            if (corelatie>ps)
            {
                ///ADAUGAREA CORELATIILOR>PS IN VECTORUL DE STRUCTURA
                (*nr)++;
                aux=(detectii*)realloc((*D),(*nr)*sizeof(detectii));
                if (aux==NULL)
                {
                    free(aux);
                    return 0;
                }
                else
                {
                    (*D)=aux;
                    (*D)[(*nr)-1].cor=corelatie;
                    (*D)[(*nr)-1].coordonate[0]=x;
                    (*D)[(*nr)-1].coordonate[1]=y;
                    (*D)[(*nr)-1].cif=NRX;
                    (*D)[(*nr)-1].ok=0;
                }
            }
            y=y+w;
            deviatie_fereastra=0;
            f_medie=0;
        }
        x=x+h;
        y=0;

    }
    for (i=0; i<h; i++)
        free(sablon[i]);
    free(sablon);
    for (i=0; i<H; i++)
        free(fereastra[i]);
    free(fereastra);

}

void contur_colorare(char *nume_imagine,char* nume_imagine_modificata,int x,int y,int w,int h,int W,int H,int NRX,contur C[])
{
    if (NRX==0)
    {
        C[NRX].culoare[0]=0;
        C[NRX].culoare[1]=0;
        C[NRX].culoare[2]=255;
    }
    if (NRX==1)
    {
        C[NRX].culoare[0]=0;
        C[NRX].culoare[1]=255;
        C[NRX].culoare[2]=255;
    }
    if (NRX==2)
    {
        C[NRX].culoare[0]=0;
        C[NRX].culoare[1]=255;
        C[NRX].culoare[2]=0;
    }
    if (NRX==3)
    {
        C[NRX].culoare[0]=255;
        C[NRX].culoare[1]=255;
        C[NRX].culoare[2]=0;
    }
    if (NRX==4)
    {
        C[NRX].culoare[0]=255;
        C[NRX].culoare[1]=0;
        C[NRX].culoare[2]=255;
    }
    if (NRX==5)
    {
        C[NRX].culoare[0]=255;
        C[NRX].culoare[1]=0;
        C[NRX].culoare[2]=0;
    }
    if (NRX==6)
    {
        C[NRX].culoare[0]=192;
        C[NRX].culoare[1]=192;
        C[NRX].culoare[2]=192;
    }
    if (NRX==7)
    {
        C[NRX].culoare[0]=0;
        C[NRX].culoare[1]=140;
        C[NRX].culoare[2]=255;
    }
    if (NRX==8)
    {
        C[NRX].culoare[0]=128;
        C[NRX].culoare[1]=0;
        C[NRX].culoare[2]=128;
    }
    if (NRX==9)
    {
        C[NRX].culoare[0]=0;
        C[NRX].culoare[1]=0;
        C[NRX].culoare[2]=128;
    }
    int i,j,k;
    pixel **imagine_modificata=img_matrice(nume_imagine_modificata,W,H);
    for (i=x; i<h+x; i++)
    {
        if (i>=H)
            break;
        j=y;
        while (j<w+y)
        {
            if (j>=W || i>=H)
                break;
            for (k=0; k<3; k++)
                imagine_modificata[i][j].octeti[k]=C[NRX].culoare[k];
            if (j==w+y-1)
            {
                j++;
                break;
            }
            if (i==x || i==h+x-1)
                j++;
            else
                j=w+y-1;

        }
    }
    unsigned char *header=Header_im(nume_imagine),a=0;
    i=0;
    FILE *out=fopen(nume_imagine_modificata,"wb");
    while (i<54)
    {
        fwrite(&header[i],sizeof(unsigned char),1,out);
        i++;
    }
    i=0;

    int padding=0;
    if (W%4!=0)
        padding=4-(3*W)%4;
    else padding=0;
    for (i=0; i<H; i++)
    {
        for (j=0; j<W; j++)
            fwrite(&imagine_modificata[i][j],sizeof(pixel),1,out);
        for (j=0; j<padding; j++)
            fwrite(&a,sizeof(unsigned char),1,out);
    }
    fclose(out);
    for (i=0; i<H; i++)
        free(imagine_modificata[i]);
    free(imagine_modificata);
    free(header);
}

void eliminare_non_maxime(detectii **D,int nr)
{
    int i,j,k,p,aux;
    double suprapunere=0,arie=0,lungime=0,latime=0;
    for (i=0; i<nr; i++)
        for (j=i+1; j<nr-1; j++)
            if ((*D)[i].cif!=(*D)[j].cif)
            {
                if ((*D)[i].coordonate[0]>(*D)[j].coordonate[0])
                {
                    k=(*D)[i].coordonate[0];
                    p=(*D)[j].coordonate[0];
                    if (k-p<=15)
                    {
                        lungime=k-p;
                        if ((*D)[i].coordonate[1]>(*D)[j].coordonate[1])
                        {
                            k=(*D)[i].coordonate[1];
                            p=(*D)[j].coordonate[1];
                            if (k-p<=11)
                            {
                                latime=k-p;
                            }
                        }
                        else
                        {
                            k=(*D)[j].coordonate[1];
                            p=(*D)[i].coordonate[1];
                            if (k-p<=11)
                                latime=k-p;
                        }
                    }
                }
                else
                {
                    k=(*D)[j].coordonate[0];
                    p=(*D)[i].coordonate[0];
                    if (k-p<=15)
                    {
                        lungime=k-p;
                        if ((*D)[i].coordonate[1]>(*D)[j].coordonate[1])
                        {
                            k=(*D)[i].coordonate[1];
                            p=(*D)[j].coordonate[1];
                            if (k-p<=11)
                            {
                                latime=k-p;
                            }
                        }
                        else
                        {
                            k=(*D)[j].coordonate[1];
                            p=(*D)[i].coordonate[1];
                            if (k-p<=11)
                                latime=k-p;
                        }
                    }


                }
                arie=lungime*latime;
                suprapunere=arie/(2*11*15-arie);
                if (suprapunere>0.2)
                {
                  (*D)[j].ok=1;
                }
            }
}

void sortare (detectii **D, int nr)
{
    qsort((*D),nr,sizeof(detectii),cmpDescrescator);
}

void coloarev2(char *nume_imagine,char *nume_imagine_modificata,detectii *D,int nr,int w,int h,int W,int H,contur C[])
{
    int i=0;
    eliminare_non_maxime(&D,nr);
    for (i=0; i<nr; i++)
    {
        contur_colorare(nume_imagine,nume_imagine_modificata,D[i].coordonate[0],D[i].coordonate[1],w,h,W,H,D[i].cif,C);
    }
}

int main()
{
    int Wp,Hp;
    char nume_imagine_initiala[30];
    char nume_imagine_criptata[40];
    char nume_imagine_decriptata[40];
    char nume_cheie_secreta[]="cheie.txt";
    printf("Introduceti numele imaginii pe care doriti sa o criptati si decriptati: \n");
    printf("EX: proiect_pp.bmp \n");
    scanf("%s",nume_imagine_initiala);
    printf("Introduceti numele pe care ati dori sa o aiba imaginea criptata: \n");
    printf("EX: criptata.bmp \n");
    scanf("%s",nume_imagine_criptata);
    printf("Introduceti numele pe care ati dori sa o aiba imaginea decriptata: \n");
    printf("EX: decriptata.bmp \n");
    scanf("%s",nume_imagine_decriptata);
    incarcare_im(nume_imagine_initiala);
    dimensiuni_imagine(nume_imagine_initiala,&Wp,&Hp);
    unsigned char *pixeli_criptati=alg_criptare(nume_imagine_initiala,nume_cheie_secreta,nume_imagine_criptata,Wp,Hp);
    creare_fisier(nume_imagine_initiala,nume_imagine_criptata,Wp,Hp,pixeli_criptati);
    unsigned char *pixeli_decriptati=alg_decriptare(nume_imagine_initiala,nume_cheie_secreta,nume_imagine_criptata,Wp,Hp,pixeli_criptati);
    creare_fisier(nume_imagine_initiala, nume_imagine_decriptata,Wp,Hp,pixeli_decriptati);
    printf("\nValorile testului chi-patrat pentru imaginea initiala %s : ",nume_imagine_initiala);
    hi_2(nume_imagine_initiala,Wp,Hp);
    printf("\nValorile testului chi_patrat pentru imaginea criptata %s : ",nume_imagine_criptata);
    hi_2(nume_imagine_criptata,Wp,Hp);
    free(pixeli_criptati);
    free(pixeli_decriptati);
    printf("PROCESUL DE CRIPTARE SI DECRIPTARE FINALIZAT!\n");



    printf("PROCESUL DE TEMPLATE-MATCHING IN DESFASURARE PE IMAGINEA test.bmp...\n");
    char nume_imagine[]="test.bmp";
    char nume_imagine_modificata[]="test1.bmp";
    grayscale_image("test.bmp", "test_grayscale.bmp");
    grayscale_image("cifra0.bmp", "cifra0_grayscale.bmp");
    grayscale_image("cifra1.bmp", "cifra1_grayscale.bmp");
    grayscale_image("cifra2.bmp", "cifra2_grayscale.bmp");
    grayscale_image("cifra3.bmp", "cifra3_grayscale.bmp");
    grayscale_image("cifra4.bmp", "cifra4_grayscale.bmp");
    grayscale_image("cifra5.bmp", "cifra5_grayscale.bmp");
    grayscale_image("cifra6.bmp", "cifra6_grayscale.bmp");
    grayscale_image("cifra7.bmp", "cifra7_grayscale.bmp");
    grayscale_image("cifra8.bmp", "cifra8_grayscale.bmp");
    grayscale_image("cifra9.bmp", "cifra9_grayscale.bmp");
    int W,H,w,h,NRX=0,i,j,nr;
    float ps=0.5;
    detectii *D=NULL;
    contur C[9];
    dimensiuni_imagine("test_grayscale.bmp",&W,&H);
    dimensiuni_imagine("cifra0_grayscale.bmp",&w,&h);
    nr=0;
    template_matching("test_grayscale.bmp","cifra0_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra1_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra2_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra3_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra4_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra5_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra6_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra7_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra8_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    NRX++;
    template_matching("test_grayscale.bmp","cifra9_grayscale.bmp",W,H,w,h,ps,NRX,C,&D,&nr);
    sortare(&D,nr);
    copie_imagine(nume_imagine,nume_imagine_modificata,W,H);
    coloarev2(nume_imagine,nume_imagine_modificata,D,nr,w,h,W,H,C);
    printf("PROCESUL DE TEMPLATE-MATCHING COMPLET!");
    return 0;

}
