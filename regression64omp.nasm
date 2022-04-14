; ---------------------------------------------------------
; Regression con istruzioni AVX a 64 bit
; ---------------------------------------------------------
; F. Angiulli
; 23/11/2017
;

;
; Software necessario per l'esecuzione:
;
;     NASM (www.nasm.us)
;     GCC (gcc.gnu.org)
;
; entrambi sono disponibili come pacchetti software 
; installabili mediante il packaging tool del sistema 
; operativo; per esempio, su Ubuntu, mediante i comandi:
;
;     sudo apt-get install nasm
;     sudo apt-get install gcc
;
; potrebbe essere necessario installare le seguenti librerie:
;
;     sudo apt-get install lib32gcc-4.8-dev (o altra versione)
;     sudo apt-get install libc6-dev-i386
;
; Per generare file oggetto:
;
;     nasm -f elf64 regression64.nasm
;

%include "sseutils64.nasm"

section .data			; Sezione contenente dati inizializzati
	u	dq	1.0
	
	
	dim	equ	8
	
	dim2 equ	4
	
	p	equ	4
	unroll equ	4
	

section .bss			; Sezione contenente dati non inizializzati

alignb 32
eta		resq	1

section .text			; Sezione contenente il codice macchina

; ----------------------------------------------------------
; macro per l'allocazione dinamica della memoria
;
;	getmem	<size>,<elements>
;
; alloca un'area di memoria di <size>*<elements> bytes
; (allineata a 16 bytes) e restituisce in EAX
; l'indirizzo del primo bytes del blocco allocato
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)
;
;	fremem	<address>
;
; dealloca l'area di memoria che ha inizio dall'indirizzo
; <address> precedentemente allocata con getmem
; (funziona mediante chiamata a funzione C, per cui
; altri registri potrebbero essere modificati)

extern get_block
extern free_block

%macro	getmem	2
	mov	rdi, %1
	mov	rsi, %2
	call	get_block
%endmacro

%macro	fremem	1
	mov	rdi, %1
	call	free_block
%endmacro


global popola
;--------------------------------------------------------------------------------------------------------------------------------
;void scriviColonne( int tuple,int coldest, int h,  int n, params* input, int* j){
;
;	int oss, i,k;
;	type num;
;	int n = input -> n;
;	int col = input-> t;
;	for(i=0;i<tuple;i++){
;			for(oss=0;oss<n;oss++){  ;J
;				num=1;
;				for(k=0;k<h;k++){
;					num=num* input->x[oss+(J[i*h+k]-1)*n];
;				}
;				input->xast[oss+coldest*n]=num;
;			}
;			coldest++; //andiamo a scrivere la colonna coldest successiva di xast
;		}
;}
;--------------------------------------------------------------------------------------------------------------------------------

;--------------------------------------------------------------------------------------------------------------------------------
;Codice con loop vectorization e loop unrolling UNROLL=4
;
 ;         for(i=0;i<tuple;i++){
;
 ;           for(oss=0;oss<n;oss=oss+p*UNROLL){
;        	type[p] num1 = [1,1,1,1];
;		type[p] num2 = [1,1,1,1];
;		type[p] num3 = [1,1,1,1];
;		type[p] num4 = [1,1,1,1];
;
;             for(k=0;k<h;k++){
;               int x_i = j[i*h+k]-1;
;              num1=num1* input->x[(oss ... oss+p-1) + j_i*n];
;             num2=num2* input->x[(oss+p ... oss+2p-1) + j_i*n];
;            num3=num3* input->x[(oss+2p ... oss+3p-1) + j_i*n];
;            num4=num4* input->x[(oss+3p ... oss+4p-1) + j_i*n];
;		}
;             input->xast[(oss ... oss+p-1)+coldest*n]=num1;
;            input->xast[(oss+p ... oss+2p-1)+coldest*n]=num2;
;           input->xast[(oss+2p ... oss+3p-1)+coldest*n]=num3;
;          input->xast[(oss+3p ... oss+4p-1)+coldest*n]=num4;
;       }
;
;           coldest++; //andiamo a scrivere la colonna coldest successiva di xast
;        }
;--------------------------------------------------------------------------------------------------------------------------


popola:
		

		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		rbp				; salva il Base Pointer
		mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
		pushaq						; salva i registri generali

		;int tuple,int coldest, int h,  int n, params* input, int* j
	
		; [R8] input->x
		; [R8+16] input->xast
		; [R8+24] input->n
		
		; RDI -> tuple
		; RSI -> coldest
		; RDX -> h
		; RCX ->  n
		; R8 -> params
		; R9 -> J
		
	;-------------------------------------------------------------------------------------------------------------
		
		
		MOV 	R10,   [R8]		; x
		MOV	        R11, [R8+16]		; xast
		
		; R8 = n/16 RBX= indice che va d a0 a R8
		MOV R8, RCX	; copio n
		SHR	R8, 4 ; n/16  numero di iterazioni del ciclo in modalità vettoriale
		

		SHL	RCX,3	; RCX = n*dim1
		SHL	RDI,2	; RDI = tuple*dim2
				
				
		;-------- i - oss - k ->  R12, R13, R14------------
		
			XOR	R12,R12	; i = 0
	
	fori:		XOR	R13,R13	; oss = 0
	
			XOR EBX, EBX	; indice
			
			CMP R8, RBX	; n/16 =0 ?
			JE	ossr		; vado nel resto
			
	foross:	
			XOR R14,R14	; k = 0
	
			VBROADCASTSD YMM0, [u]		; YMM0= [1.0, 1.0, 1.0, 1.0]	num1
			VBROADCASTSD YMM2, [u]		; YMM2= [1.0, 1.0, 1.0, 1.0]	num2
			VBROADCASTSD YMM4, [u]		; YMM4= [1.0, 1.0, 1.0, 1.0]	num3
			VBROADCASTSD YMM6, [u]		; YMM6= [1.0, 1.0, 1.0, 1.0]	num4
			
	fork:		
			;-------------------------J[i*h+k]--------------------------------------------------------------------------------------
			
			MOV		RAX, R12		; RAX =i*dim2	
			IMUL	RAX, RDX		; RAX = i*dim2*h
		
			ADD		RAX, R14		; RAX= i*h*dim2 + k*dim2
		
			MOV		EAX, [R9+RAX]	; j_i 
		
		
			SUB		EAX,1		; j_i - 1
			
			;---------------------x[oss+(j_i - 1)*n]------------------------------------------------------------------------
		
		
			IMUL	RAX, RCX		; (j_i - 1)* n
			ADD		RAX, R13		; (j_i - 1)* n + oss
			
			
			VMOVUPD	YMM1, [R10 + RAX]	; YMM1 = input->x[(oss ... oss+p-1) + j_i*n];
			VMULPD	YMM0,YMM1			; num1=num1* input->x[(oss ... oss+p-1) + j_i*n];
			
			VMOVUPD	YMM3, [R10 + RAX + 32]	; YMM3 = input->x[(oss ... oss+p-1) + j_i*n];
			VMULPD	YMM2,YMM3			; num2=num2* input->x[(oss ... oss+p-1) + j_i*n];
			
			VMOVUPD	YMM5, [R10 + RAX + 64]	; YMM5 = input->x[(oss ... oss+p-1) + j_i*n];
			VMULPD	YMM4,YMM5			; num3=num3* input->x[(oss ... oss+p-1) + j_i*n];
			
			VMOVUPD	YMM7, [R10 + RAX + 96]	; YMM7 = input->x[(oss ... oss+p-1) + j_i*n];
			VMULPD	YMM6,YMM7			; num4=num4* input->x[(oss ... oss+p-1) + j_i*n];
			
			
			ADD 	R14, dim2		;k++ 	k=k+dim2
			IMUL	R15, RDX, dim2	;R14= h*dim2
			CMP		R14, R15			;k<h?
			JL		fork
			
			
			;------------------input->xast[(oss ... oss+p-1)+coldest*n]=num----------------------------------------------
		
			
			
			MOV		RAX, RSI		; RBX = coldest
			IMUL	RAX, RCX		; EBX = coldest*n*dim
			ADD		RAX, R13		; R14 = oss*dim + coldest*n*dim
		
			VMOVUPD	[R11 + RAX], YMM0		; xast[(oss ... oss+p-1)+coldest*n]=num1;
			VMOVUPD	[R11 + RAX + 32], YMM2	; xast[(oss+p ... oss+2p-1)+coldest*n]=num2;
			VMOVUPD	[R11 + RAX + 64], YMM4	; xast[(oss+2p ... oss+3p-1)+coldest*n]=num3;
			VMOVUPD	[R11 + RAX + 96], YMM6	; xast[(oss+3p ... oss+4p-1)+coldest*n]=num4;
			
			ADD 	R13, dim*p*unroll	; oss++
			
			ADD 	RBX, 1
			
			CMP 	RBX, R8		
			JL 		foross
			
			CMP		R13, RCX		; oss=n?
			JE		incri
			
			;-----------------------------------------------------for resto---------------------------------------------------------------------	
			
	ossr:		XOR R14,R14	; k = 0
			VBROADCASTSD YMM0, [u]		; YMM1= [1.0, 1.0, 1.0, 1.0]
	
	forkr:	
			;-------------------------J[i*h+k]--------------------------------------------------------------------------------------
			
			MOV		RAX, R12		; RAX =i*dim2	
			IMUL	RAX, RDX		; RAX = i*dim2*h
		
			ADD		RAX, R14		; RAX= i*h*dim2 + k*dim2
		
			MOV		EAX, [R9+RAX]	; j_i 
		
		
			SUB		EAX,1		; j_i - 1
			
			;---------------------x[oss+(j_i - 1)*n]------------------------------------------------------------------------
		
		
			IMUL	RAX, RCX		; (j_i - 1)* n
			ADD		RAX, R13		; (j_i - 1)* n + oss
			
			
			VMOVSD	XMM1, [R10 + RAX]	; XMM1 = input->x[(oss+ x_i*n];
			VMULSD	XMM0,XMM1			; num1=num1* input->x[(oss  + x_i*n];
			
		
			ADD 	R14, dim2		;k++ 	k=k+dim2
			IMUL	R15, RDX, dim2	;R14= h*dim2
			CMP		R14, R15			;k<h?
			JL		forkr
	
			
			;------------------------------input->xast[oss+coldest*n]=num;----------------------------
		
			MOV		RAX, RSI		; RBX = coldest
			IMUL	RAX, RCX		; EBX = coldest*n*dim
			ADD		RAX, R13		; R14 = oss*dim + coldest*n*dim
			
			
			VMOVSD	[R11 + RAX], XMM0
			
			ADD 	R13, dim		; oss++
			CMP		R13, RCX		; oss<n?
			JL 		ossr
			
		
	incri:		ADD		RSI,1		; coldest++
			ADD 	R12, dim2	; i++
			CMP		R12, RDI		; i<tuple?
			JL		fori
			
			

		
		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------
		
		popaq						; ripristina i registri generali
		mov	rsp, rbp				; ripristina lo Stack Pointer
		pop		rbp					; ripristina il Base Pointer
		ret							; torna alla funzione C chiamante



;------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


		
global sommaVettori

sommaVettori:

				push		rbp				; salva il Base Pointer
				mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
				pushaq						; salva i registri generali

				;[RDI]  daSottrarre
				;[RSI]   temp
				;RDX    size
				
				MOV		RAX, 0			;i
				MOV		R8, RDX			;copia per ciclo fors
				SHR			RDX, 4		;iterazioni possibili a gruppi di 4
				MOV		RCX, 0			;calcola lo sfasamento di 32 byte nel for e di 8 byte nel fors
			for:	CMP		RAX, RDX
				JNL			fine
				MOV		RCX, RAX
				SHL			RCX, 7
				VMOVUPD 	YMM0, [RDI+RCX]
				VMOVUPD 	YMM1, [RDI+RCX+32]
				VMOVUPD 	YMM2, [RDI+RCX+64]
				VMOVUPD 	YMM3, [RDI+RCX+96]
				VMOVUPD	YMM4, [RSI+RCX]
				VMOVUPD	YMM5, [RSI+RCX+32]
				VMOVUPD	YMM6, [RSI+RCX+64]
				VMOVUPD	YMM7, [RSI+RCX+96]
				VADDPD		YMM0, YMM4
				VADDPD		YMM1, YMM5
				VADDPD		YMM2, YMM6
				VADDPD		YMM3, YMM7
				VMOVUPD	[RDI+RCX], YMM0
				VMOVUPD	[RDI+RCX+32], YMM1
				VMOVUPD	[RDI+RCX+64], YMM2
				VMOVUPD	[RDI+RCX+96], YMM3
				ADD		RAX, 1
				JMP			for
			fine: SHL			RDX, 4		;calcola da quale indice i ricominciare
			fors:	CMP		RDX, R8
				JNL			return
				MOV		RCX, RDX
				SHL			RCX, 3
				VMOVSD 	XMM0, [RDI+RCX]
				VMOVSD		XMM1, [RSI+RCX]
				VADDSD		XMM0, XMM1
				VMOVSD		[RDI+RCX], XMM0
				ADD		RDX, 1
				JMP			fors
			return:			popaq						; ripristina i registri generali
							mov	rsp, rbp			; ripristina lo Stack Pointer
							pop		rbp					; ripristina il Base Pointer
							ret							; torna alla funzione C chiamante
							
							
global moltiplicazionePerScalare

moltiplicazionePerScalare:

						push		rbp				; salva il Base Pointer
						mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
						pushaq						; salva i registri generali
						
						;[RDI]  vett_ris
						;[RSI]  x
						;RDX   i*size
						;RCX   size
						;XMM0  scalare
						
						VBROADCASTSD	YMM2, XMM0
						MOV	RAX, 0
						MOV	R8, RCX
						SHR		RCX, 4
						SHL		RDX, 3
					for1 :	CMP	RAX, RCX
						JNL		fine1
						MOV	R9, RAX
						SHL		R9, 7
						ADD	RDX, R9
						VMOVUPD	YMM1, [RSI+RDX]
						VMOVUPD	YMM3, [RSI+RDX+32]
						VMOVUPD	YMM4, [RSI+RDX+64]
						VMOVUPD	YMM5, [RSI+RDX+96]
						VMULPD		YMM1, YMM2
						VMULPD		YMM3, YMM2
						VMULPD		YMM4, YMM2
						VMULPD		YMM5, YMM2
						VMOVUPD	[RDI+R9], YMM1
						VMOVUPD	[RDI+R9+32], YMM3
						VMOVUPD	[RDI+R9+64], YMM4
						VMOVUPD	[RDI+R9+96], YMM5
						SUB		RDX, R9
						ADD	RAX, 1
						JMP		for1
					fine1:	SHL		RCX, 4
					fors1:	CMP	RCX, R8
						JNL		return1
						MOV	R9, RCX
						SHL		R9, 3
						ADD	RDX, R9
						VMOVSD	XMM1, [RSI+RDX]
						VMULSD	XMM1, XMM2
						VMOVSD	[RDI+R9], XMM1
						SUB		RDX, R9
						ADD	RCX, 1
						JMP		fors1
					return1:	popaq						; ripristina i registri generali
						mov	rsp, rbp			; ripristina lo Stack Pointer
						pop		rbp					; ripristina il Base Pointer
						ret							; torna alla funzione C chiamante
						
						
global quadratoVettore

quadratoVettore:

						push		rbp				; salva il Base Pointer
						mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
						pushaq						; salva i registri generali
						
						;[RDI] g
						;RSI   size*riga
						;[RDX] G
						;RCX   size
						
						MOV	RAX, 0
						MOV	R8, RCX
						SHR		RCX, 4
						SHL		RSI, 3
					for2:CMP	RAX, RCX
						JNL		fine2
						MOV	R9, RAX
						SHL		R9, 7
						VMOVUPD	YMM0, [RDI+R9]
						VMOVUPD	YMM1, [RDI+R9+32]
						VMOVUPD	YMM2, [RDI+R9+64]
						VMOVUPD	YMM3, [RDI+R9+96]
						ADD	R9, RSI
						VFMADD213PD	YMM0,YMM0,[RDX+R9]
						VFMADD213PD	YMM1,YMM1,[RDX+R9+32]
						VFMADD213PD	YMM2,YMM2,[RDX+R9+64]
						VFMADD213PD	YMM3,YMM3,[RDX+R9+96]
						VMOVUPD		[RDX+R9], YMM0
						VMOVUPD		[RDX+R9+32], YMM1
						VMOVUPD		[RDX+R9+64], YMM2
						VMOVUPD		[RDX+R9+96], YMM3
						ADD	RAX, 1
						JMP		for2
					fine2:SHL	RCX, 4
					fors2:	CMP	RCX, R8
						JNL		return2
						MOV	R9, RCX
						SHL		R9, 3
						VMOVSD	XMM0, [RDI+R9]
						ADD	R9, RSI
						VFMADD213SD	XMM0,XMM0,[RDX+R9]
						VMOVSD		[RDX+R9], XMM0
						ADD	RCX, 1
						JMP		fors2
					return2:		popaq						; ripristina i registri generali
						mov	rsp, rbp			; ripristina lo Stack Pointer
						pop		rbp					; ripristina il Base Pointer
						ret							; torna alla funzione C chiamante
						
global faiSommatoria

faiSommatoria:

						push		rbp				; salva il Base Pointer
						mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
						pushaq						; salva i registri generali
						
						;[RDI]  daSottrarre
						;RSI	    size
						;[RDX] g
						;[RCX] G
						;R8	  riga*size
						;XMM0  epsilon
						VBROADCASTSD	YMM2, XMM0
						MOV	RAX, 0
						MOV	R9, RSI
						SHR		RSI, 2
						SHL		R8, 3
					for3:CMP	RAX,RSI
						JNL		fine3
						MOV	R10, RAX
						SHL		R10, 5
						ADD	R10, R8
						VSQRTPD	YMM1, [RCX+R10]
						VADDPD		YMM1, YMM2
						SUB		R10, R8
						VMOVUPD	YMM3, [RDX+R10]
						VDIVPD		YMM3, YMM1
						VMOVUPD	YMM4, [RDI+R10]
						VADDPD		YMM4, YMM3
						VMOVUPD	[RDI+R10], YMM4
						ADD	RAX, 1
						JMP		for3
					fine3: SHL	RSI, 2
					fors3: CMP	RSI, R9
						JNL		return3
						MOV	R10, RSI
						SHL		R10, 3
						ADD	R10, R8
						VSQRTSD	XMM1, [RCX+R10]
						VADDSD		XMM1, XMM2
						SUB		R10, R8
						VMOVSD	XMM3, [RDX+R10]
						VDIVSD		XMM3, XMM1
						VMOVSD		XMM4, [RDI+R10]
						VADDSD		XMM4, XMM3
						VMOVSD		[RDI+R10], XMM4
						ADD	RSI, 1
						JMP		fors3
					return3: popaq						; ripristina i registri generali
						mov	rsp, rbp			; ripristina lo Stack Pointer
						pop		rbp					; ripristina il Base Pointer
						ret							; torna alla funzione C chiamante
						
global prodottoScalare

prodottoScalare:

						push		rbp				; salva il Base Pointer
						mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
						pushaq						; salva i registri generali
						
						;[RDI]   cost
						;[RSI]   theta
						;[RDX]  xast
						;RCX  j*size
						;R8    size
						
						MOV	RAX, 0
						MOV	R9,R8
						SHR		R8, 4
						SHL		RCX, 3
						VXORPD	YMM0, YMM0
					for4:CMP	RAX, R8
						JNL		fine4
						MOV	R10, RAX
						SHL		R10, 7
						VMOVUPD	YMM1, [RSI+R10]
						VMOVUPD	YMM2, [RSI+R10+32]
						VMOVUPD	YMM3, [RSI+R10+64]
						VMOVUPD	YMM4, [RSI+R10+96]
						ADD	R10, RCX
						VFMADD231PD	YMM0, YMM1, [RDX+R10]
						VFMADD231PD	YMM0, YMM2, [RDX+R10+32]
						VFMADD231PD	YMM0, YMM3, [RDX+R10+64]
						VFMADD231PD	YMM0, YMM4, [RDX+R10+96]
						ADD	RAX, 1
						JMP		for4
					fine4: SHL	R8, 4
						;VXORPD	XMM3, XMM3
						VHADDPD	YMM0, YMM0
						VPERMPD	YMM0, YMM0, 11011000b
						VHADDPD	YMM0, YMM0
					fors4: CMP	R8, R9
						JNL		return4
						MOV	R10, R8
						SHL		R10, 3
						VMOVSD	XMM1, [RSI+R10]
						ADD	R10, RCX
						VFMADD231SD	XMM0, XMM1, [RDX+R10]
						ADD	R8, 1
						JMP		fors4
					return4: VMOVSD	[RDI], XMM0
						popaq						; ripristina i registri generali
						mov	rsp, rbp			; ripristina lo Stack Pointer
						pop		rbp					; ripristina il Base Pointer
						ret							; torna alla funzione C chiamante
						
						
global calcolaTheta

calcolaTheta:

						push		rbp				; salva il Base Pointer
						mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
						pushaq						; salva i registri generali
						
						;[RDI]    daSottrarre
						;[RSI]    theta
						;XMM0  eta/v
						;RDX     size

						VBROADCASTSD	YMM2, XMM0
						MOV 	RAX, 0
						MOV	RCX, RDX
						SHR		RCX, 4
					for5: CMP	RAX, RCX
						JNL		fine5
						SHL		RAX, 7
						VMOVUPD	YMM1, [RDI+RAX]
						VMOVUPD	YMM3, [RDI+RAX+32]
						VMOVUPD	YMM4, [RDI+RAX+64]
						VMOVUPD	YMM5, [RDI+RAX+96]
						VFNMADD213PD	YMM1, YMM2, [RSI+RAX]
						VFNMADD213PD	YMM3, YMM2, [RSI+RAX+32]
						VFNMADD213PD	YMM4, YMM2, [RSI+RAX+64]
						VFNMADD213PD	YMM5, YMM2, [RSI+RAX+96]
						VMOVUPD	[RSI+RAX], YMM1
						VMOVUPD	[RSI+RAX+32], YMM3
						VMOVUPD	[RSI+RAX+64], YMM4
						VMOVUPD	[RSI+RAX+96], YMM5
						SHR		RAX, 7
						ADD	RAX, 1
						JMP		for5
					fine5: SHL	RCX, 4
					fors5: CMP	RCX, RDX
						JNL		return5
						SHL		RCX, 3
						VMOVSD	XMM1, [RDI+RCX]
						VFNMADD213SD	XMM1, XMM2, [RSI+RCX]
						VMOVSD	[RSI+RCX], XMM1
						SHR		RCX, 3
						ADD	RCX, 1
						JMP		fors5
					return5:	popaq						; ripristina i registri generali
						mov	rsp, rbp			; ripristina lo Stack Pointer
						pop		rbp					; ripristina il Base Pointer
						ret							; torna alla funzione C chiamante
						
						
global azzeramentoDaSottrarre

azzeramentoDaSottrarre:

						push		rbp				; salva il Base Pointer
						mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
						pushaq						; salva i registri generali
						
						;[RDI]  daSottrarre
						;RSI	    size
						
						VXORPD	YMM0, YMM0
						MOV	RAX, 0
						MOV	RCX, RSI
						SHR		RCX, 4
					for6:CMP	RAX, RCX
						JNL		fine6
						SHL		RAX, 7
						VMOVUPD	[RDI+RAX], YMM0
						VMOVUPD	[RDI+RAX+32], YMM0
						VMOVUPD	[RDI+RAX+64], YMM0
						VMOVUPD	[RDI+RAX+96], YMM0
						SHR		RAX, 7
						ADD	RAX, 1
						JMP		for6
					fine6: SHL	RCX, 4
					fors6: CMP	RCX, RSI
						JNL		return6
						SHL		RCX, 3
						VMOVSD	[RDI+RCX], XMM0
						SHR		RCX, 3
						ADD	RCX, 1
						JMP		fors6
					return6:		popaq						; ripristina i registri generali
						mov	rsp, rbp			; ripristina lo Stack Pointer
						pop		rbp					; ripristina il Base Pointer
						ret							; torna alla funzione C chiamante