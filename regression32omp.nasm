;---------------------------------------------------------
; Regressione con istruzioni SSE a 32 bit
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
;     nasm -f elf32 regression32.nasm 
;
%include "sseutils32.nasm"

section .data			; Sezione contenente dati inizializzati
	dim 	equ 4
	p 		equ 4
	unroll	equ 4
	align 16
	uni		dd	1.0, 1.0, 1.0, 1.0


section .bss			; Sezione contenente dati non inizializzati
	alignb 16
	eta		resd		1

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
	mov	eax, %1
	push	eax
	mov	eax, %2
	push	eax
	call	get_block
	add	esp, 8
%endmacro

%macro	fremem	1
	push	%1
	call	free_block
	add	esp, 4
%endmacro

; ------------------------------------------------------------
; Funzioni
; ------------------------------------------------------------

;-------------------------------------------------------------------------
; Funzione popola
; ------------------------------------------------------------------------------

global popola1


;for(i=0;i<tuple;i++){
;
;	for(oss=0;oss<n;oss=oss+p*unroll){
;		type[p] xmm0 = [1,1,1,1]
;		type[p] xmm2 = [1,1,1,1]
;		type[p] xmm4 = [1,1,1,1]
;		type[p] xmm6 = [1,1,1,1]
;		for(k=0;k<h;k++){
;			int x_i = j[i*h+k]
;			num1=num1* ->x[(oss ... oss+p-1) + x_i*n]
;			num1=num1* ->x[(oss+p ... oss+2p-1) + x_i*n]
;			num1=num1* ->x[(oss+2p ... oss+3p-1) + x_i*n]
;			num1=num1* ->x[(oss+3p ... oss+4p-1) + x_i*n]
;		}
;		input->xast[(oss ... oss+p-1)+coldest*n]=num1
;		input->xast[(oss ... oss+p-1)+coldest*n]=num1
;		input->xast[(oss ... oss+p-1)+coldest*n]=num1
;		input->xast[(oss ... oss+p-1)+coldest*n]=num1
;	}
;	for(oss; oss<n; oss++) RESTO
;	coldest++; //andiamo a scrivere la colonna coldest successiva di xast
;}



tuple 	equ 	8
n		equ		12
coldest	equ		16
h		equ		20
d		equ		24
j		equ		28
col		equ		32
x		equ		36
xast 	equ		40
popola1:
		; ------------------------------------------------------------
		; Sequenza di ingresso nella funzione
		; ------------------------------------------------------------
		push		ebp							; salva il Base Pointer
		mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
		push		ebx							; salva i registri da preservare
		push		esi
		push		edi
		
		
		
	
		
		; [EcX] input->x
		; [EcX+4] input->y
		; [EcX+8] input->xast
		; [EcX+12] input->n
		; [EAX+16] input->d
		; [EcX+20] input->k
		; [EcX+24] input->degree
		; [EcX+28] input->eta

		mov			ecx, [ebp+n]
		
		
		xor			esi,esi						;esi = i = 0
fori:	xor			edi, edi					;edi = oss = 0---> ricordo oss lo devo incrementare di p
foross:	xor			edx, edx					;edx = k = 0
		MOVAPS		XMM0, [uni]
		MOVAPS		XMM2, [uni]
		MOVAPS		XMM4, [uni]
		MOVAPS		XMM6, [uni]
fork:	mov			eax, [ebp+j]				;eax = j
		mov			ebx, [ebp+h]				;ebx = h
		imul		ebx, esi					;ebx = i*h
		add			ebx, edx					;ebx = i*h+k
; -------------------------------------------------  num1=num1* ->x[(oss ... oss+p-1) + x_i*n]
		mov			eax, [eax+ebx*dim]			;X_I -> j[i*h*dim + k*dim]
		dec			eax 	
		imul		eax, [ebp+n]				;x_i * n	
		add			eax, edi					; x_i * n + oss
		mov			ebx, [ebp+x]				;ebx = input->x

		MOVUPS		XMM1,[ebx + eax*dim] 		; XMM1 = input->x[(oss ... oss+p-1) + x_i*n];
		MULPS		XMM0, XMM1					; num1=num1* input->x[(oss ... oss+p-1) + x_i*n];

		MOVUPS		XMM3,[ebx + 16 + eax*dim] 	; XMM1 = input->x[(oss+p ... oss+2p-1) + x_i*n];
		MULPS		XMM2, XMM3					; num1=num1* input->x[(oss+p ... oss+2p-1) + x_i*n];
		
		MOVUPS		XMM5,[ebx + 32 + eax*dim] 	; XMM1 = input->x[(oss+2p ... oss+3p-1) + x_i*n];
		MULPS		XMM4, XMM5					; num1=num1* input->x[(oss+2p... oss+3p-1) + x_i*n];

		MOVUPS		XMM7,[ebx + 48 + eax*dim] 	; XMM1 = input->x[(oss+3p ... oss+4p-1) + x_i*n];
		MULPS		XMM6, XMM7					; num1=num1* input->x[(oss+3p ... oss+4p-1) + x_i*n];

		inc			edx							; k++ 
		cmp			edx, [ebp+h]
		jl			fork
;------------------------------------------------  input->xast[(oss ... oss+p-1)+coldest*n]=num1

		mov			ebx, [ebp+coldest]			;ebx=coldest
		add			ebx, esi					;coldest giusta (perché non faccio coldest++)
		imul		ebx, [ebp+n]				;ebx=coldest*n
		add			ebx, edi      				;ebx=coldest*n +oss
		mov			eax, [ebp+xast]				; eax=x_ast
		MOVUPS 		[eax+ebx*dim], XMM0			;x_ast[oss*dim + coldest*n*dim] = xmm0
		MOVUPS 		[eax+16+ebx*dim], XMM2		;x_ast[oss*dim + coldest*n*dim] = xmm2
		MOVUPS 		[eax+32+ebx*dim], XMM4		;x_ast[oss*dim + coldest*n*dim] = xmm4	
		MOVUPS 		[eax+48+ebx*dim], XMM6		;x_ast[oss*dim + coldest*n*dim] = xmm6
		add			edi, 16
		mov			eax, ecx					;eax=n
		sub 		eax, edi					;n-oss
		cmp			eax, 16
		jge			foross	
		cmp			edi, ecx					;oss==n?
		je			nr


;----------------------------------------------------for resto
ossr:	xor			edx, edx					;edx = k = 0
		MOVAPS		XMM0, [uni]
forkr:	mov			eax, [ebp+j]				;eax = j
		mov			ebx, [ebp+h]				;ebx = h
		imul		ebx, esi					;ebx = i*h
		add			ebx, edx					;ebx = i*h+k
; -------------------------------------------------  num1=num1* ->x[(oss ... oss+p-1) + x_i*n]
		mov			eax, [eax+ebx*dim]			;X_I -> j[i*h*dim + k*dim]
		dec			eax
		imul		eax, [ebp+n]				;x_i * n	
		add			eax, edi					; x_i * n + oss
		mov			ebx, [ebp+x]				;ebx = input->x
		MOVSS		XMM1,[ebx + eax*dim] 		; XMM1 = input->x[(oss+ x_i*n];
		MULSS		XMM0, XMM1					; num1=num1* input->x[(oss  + x_i*n];
		inc			edx							; k++ 
		cmp			edx, [ebp+h]
		jl			forkr
;------------------------------------------------  input->xast[(oss ... oss+p-1)+coldest*n]=num1

		mov			ebx, [ebp+coldest]			;ebx=coldest
		add			ebx, esi					;coldest giusta (perché non faccio coldest++)
		imul		ebx, [ebp+n]				;ebx=coldest*n
		add			ebx, edi					;ebx=coldest*n +oss
		mov			eax, [ebp+xast]				; eax=x_ast
		MOVSS 		[eax+ebx*dim], XMM0			;x_ast[oss*dim + coldest*n*dim] = xmm1
		inc 		edi							;oss++
		cmp			edi, ecx			
		jl	 		ossr
;-------------------------------------------------------fine for resto
nr:		inc			esi							;i++
		cmp			esi, [ebp+tuple]
		jl			fori
		

		
		




		; ------------------------------------------------------------
		; Sequenza di uscita dalla funzione
		; ------------------------------------------------------------

		pop	edi									; ripristina i registri da preservare
		pop	esi
		pop	ebx
		mov	esp, ebp							; ripristina lo Stack Pointer
		pop	ebp									; ripristina il Base Pointer
		ret	






global prodottoScalare

prodottoScalare:

				push		ebp							; salva il Base Pointer
				mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
				push		ebx							; salva i registri da preservare
				push		esi
				push		edi
		
				
				MOV	 EDX, [EBP+12] ;theta
				MOV 	EBX, [EBP+16] ;xast
				MOV 	EDI,  [EBP+20] ;j*size
				MOV 	ESI,  [EBP+24] ;size
				
				
				MOV 	EAX, 0	;indice i
				
				MOV 	ECX, ESI ;servirà per il for scalare
				
				SHR	 	ESI, 2		;calcolo size//4
				SHL  	EDI, 2		;calcolo j*size*dim
				XORPS	XMM0,XMM0
				
			for:	CMP	 	EAX, ESI
				JNL		 fine
				SHL		EAX, 4				;i*16
				MOVUPS	XMM1, [EDX+EAX]
				ADD	 	EAX, EDI      			;j*size+i
				MOVUPS	XMM2, [EBX+EAX]
				SUB	 	EAX, EDI
				MULPS	XMM1,XMM2
				ADDPS	XMM0, XMM1
				SHR		EAX, 4				;i/16
				ADD		EAX, 1				;i++
				JMP	 	for
				
			fine: SHL	 	ESI, 2
			
			fors: CMP 	ESI, ECX
				JNL	 	return
				SHL		ESI, 2
				MOVSS	XMM1, [EDX+ESI]
				ADD	 	ESI, EDI
				MOVSS	XMM2, [EBX+ESI]
				SUB	 	ESI, EDI
				MULSS	XMM1,XMM2
				ADDSS	XMM0, XMM1
				SHR		ESI, 2
				ADD		ESI, 1
				JMP	 	fors
				
			return:	HADDPS	XMM0, XMM0
					HADDPS	XMM0, XMM0
					MOV	EAX, [EBP+8]
					MOVSS	[EAX], XMM0
	
				
				pop	edi									; ripristina i registri da preservare
				pop	esi
				pop	ebx
				mov	esp, ebp							; ripristina lo Stack Pointer
				pop	ebp									; ripristina il Base Pointer
				ret										; torna alla funzione C chiamante
				

		
global sommaVettori

sommaVettori:

				push		ebp							; salva il Base Pointer
				mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
				push		ebx							; salva i registri da preservare
				push		esi
				push		edi
		
				MOV	EAX,[EBP+8]		;daSottrarre
				MOV	EBX,[EBP+12]		;temp
				MOV	EDX,[EBP+16]		;size
				
				SHR		EDX,4			;size//4
				
				XOR		EDI,EDI
				
				
			for1:	CMP		EDI,EDX
				JNL		fine1
				SHL		EDI,6
				
				MOVUPS	XMM0,[EAX+EDI]
				MOVUPS	XMM1,[EAX+EDI+16]
				MOVUPS	XMM2,[EAX+EDI+32]
				MOVUPS	XMM3,[EAX+EDI+48]
				
				MOVUPS	XMM4,[EBX+EDI]
				MOVUPS	XMM5,[EBX+EDI+16]
				MOVUPS	XMM6,[EBX+EDI+32]
				MOVUPS	XMM7,[EBX+EDI+48]
				
				ADDPS	XMM0,XMM4
				ADDPS	XMM1,XMM5
				ADDPS	XMM2,XMM6
				ADDPS	XMM3,XMM7
				
				MOVUPS	[EAX+EDI],XMM0
				MOVUPS	[EAX+EDI+16],XMM1
				MOVUPS	[EAX+EDI+32],XMM2
				MOVUPS	[EAX+EDI+48],XMM3
				
				SHR		EDI,6
				ADD		EDI,1
				JMP		for1
				
			fine1:	SHL		EDX,4
		
				
			fors1: CMP		EDI,EDX
				JNL		return1
				SHL		EDI,2
				MOVSS	XMM1,[EAX+EDI]
				MOVSS	XMM2,[EBX+EDI]
				ADDSS	XMM1,XMM2
				MOVSS	[EAX+EDI],XMM1
				SHR		EDI,2
				ADD		EDI,1
				JMP		fors1
				
	
				
			return1:	pop	edi									; ripristina i registri da preservare
					pop	esi
					pop	ebx
					mov	esp, ebp							; ripristina lo Stack Pointer
					pop	ebp									; ripristina il Base Pointer
					ret										; torna alla funzione C chiamante

global moltiplicazionePerScalare

moltiplicazionePerScalare:

				push		ebp							; salva il Base Pointer
				mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
				push		ebx							; salva i registri da preservare
				push		esi
				push		edi
		
				MOV	EAX,[EBP+8]		;vettris
				MOV	EBX,[EBP+12]		;Matrix x
				MOV	ECX,[EBP+16]		;i*size
				MOV	EDX,[EBP+20]		;size
				MOVSS	XMM0,[EBP+24]	;scalare
				SHUFPS	XMM0,XMM0,00000000
				
				MOV	EDI,EDX			;copia di size
				XOR		ESI,ESI
				
		
				SHR		EDX,4			;size//4
				SHL		ECX,2			;i*size*dim
				
				for2:		CMP		ESI,EDX
						JNL		fine2
						SHL		ESI,6
						
						ADD		ECX,ESI
						
						MOVUPS	XMM2, [EBX+ECX]
						MOVUPS	XMM3, [EBX+ECX+16]
						MOVUPS	XMM4, [EBX+ECX+32]
						MOVUPS	XMM5, [EBX+ECX+48]
						
						MULPS	XMM2,XMM0
						MULPS	XMM3,XMM0
						MULPS	XMM4,XMM0
						MULPS	XMM5,XMM0
						
						MOVUPS	[EAX+ESI], XMM2
						MOVUPS	[EAX+ESI+16], XMM3
						MOVUPS	[EAX+ESI+32], XMM4
						MOVUPS	[EAX+ESI+48], XMM5
						
						SUB		ECX,ESI
						
						SHR		ESI,6
						ADD		ESI,1
						JMP		for2
						
				fine2:	SHL		EDX,4
						
						
				fors2:	CMP		EDX,EDI
						JNL		return2
						SHL		EDX,2
						ADD		ECX,EDX
						MOVSS	XMM3,[EBX+ECX]
						SUB		ECX,EDX
						MULSS	XMM3,XMM0
						MOVSS	[EAX+EDX],XMM3
						SHR		EDX,2
						ADD		EDX,1
						JMP		fors2
						
				
				return2:	pop	edi									; ripristina i registri da preservare
						pop	esi
						pop	ebx
						mov	esp, ebp								; ripristina lo Stack Pointer
						pop	ebp									; ripristina il Base Pointer
						ret										; torna alla funzione C chiamante



global quadratoVettore

quadratoVettore:

				push		ebp							; salva il Base Pointer
				mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
				push		ebx							; salva i registri da preservare
				push		esi
				push		edi
		
				MOV	EAX,[EBP+8]		;vector g
				MOV	EBX,[EBP+12]		;riga*size
				MOV	ECX,[EBP+16]		;Matrix G
				MOV	EDX,[EBP+20]		;size
				
				XOR		ESI,ESI
				MOV	EDI,EDX			;copia di size
		
				SHR		EDX,4			;size//4
				SHL		EBX,2			;riga*size*dim
				
				for3:		CMP		ESI,EDX
						JNL		fine3
						SHL		ESI,6
						
						MOVUPS	XMM0,[EAX+ESI]		;g[i]
						MOVUPS	XMM1,[EAX+ESI+16]	
						MOVUPS	XMM2,[EAX+ESI+32]	
						MOVUPS	XMM3,[EAX+ESI+48]	
						
						MULPS	XMM0,XMM0			;g[i]*g[i]
						MULPS	XMM1,XMM1	
						MULPS	XMM2,XMM2	
						MULPS	XMM3,XMM3	
						
						ADD		ESI,EBX
						
						MOVUPS	XMM4, [ECX+ESI]		;G[riga*size+i]
						MOVUPS	XMM5, [ECX+ESI+16]
						MOVUPS	XMM6, [ECX+ESI+32]
						MOVUPS	XMM7, [ECX+ESI+48]
						
						ADDPS	XMM4,XMM0
						ADDPS	XMM5,XMM1
						ADDPS	XMM6,XMM2
						ADDPS	XMM7,XMM3
						
						MOVUPS	[ECX+ESI], XMM4
						MOVUPS	[ECX+ESI+16], XMM5
						MOVUPS	[ECX+ESI+32], XMM6
						MOVUPS	[ECX+ESI+48], XMM7
						
						SUB		ESI,EBX
						
						SHR		ESI,6
						ADD		ESI,1
						JMP		for3
						
				fine3:	SHL		EDX,4
						
						
				fors3:	CMP		EDX,EDI
						JNL		return3
						SHL		EDX,2
						MOVSS	XMM0,[EAX+EDX]		;g[i]
						MULSS	XMM0,XMM0			;g[i]*g[i]
						ADD		EDX,EBX
						MOVSS	XMM1, [ECX+EDX]		;G[riga*size+i]
						ADDSS	XMM1,XMM0
						MOVSS	[ECX+EDX], XMM1
						SUB		EDX,EBX
						SHR		EDX,2
						ADD		EDX,1
						JMP		fors3
						
				
				return3:	pop	edi									; ripristina i registri da preservare
						pop	esi
						pop	ebx
						mov	esp, ebp								; ripristina lo Stack Pointer
						pop	ebp									; ripristina il Base Pointer
						ret										; torna alla funzione C chiamante

global faiSommatoria

faiSommatoria:

				push		ebp							; salva il Base Pointer
				mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
				push		ebx							; salva i registri da preservare
				push		esi
				push		edi
				
				
				MOV	EAX,[EBP+8]		;vector daSottrarre
				MOV	EBX,[EBP+12]		;size
				MOV	ECX,[EBP+16]		;vector g
				MOV	EDX,[EBP+20]		;matrix G
				MOV	EDI,[EBP+24]		;riga*size
				MOVSS	XMM0,[EBP+28]	;epsilon
				SHUFPS	XMM0,XMM0,00000000
				
				
				XOR		ESI,ESI
				;MOV	EDI,EDX			;copia di size
		
				SHR		EBX,2			;size//4
				SHL		EDI,2			;riga*size*dim
				
				for4:		CMP		ESI,EBX
						JNL		fine4
						SHL		ESI,4				;h*dim
						MOVUPS	XMM1,[ECX+ESI]		;g[h]
						
						
						
						ADD		ESI,EDI
						MOVUPS	XMM2, [EDX+ESI]
						SQRTPS	XMM2,XMM2	;radice G[riga*size+h]
						
						SUB		ESI,EDI
						ADDPS	XMM2,XMM0			;radice+epsilon
						DIVPS	XMM1,XMM2
						
						MOVUPS	XMM3,[EAX+ESI]		;daSottrarre[h]
						ADDPS	XMM3,XMM1
						MOVUPS	[EAX+ESI],XMM3
						SHR		ESI,4
						ADD		ESI,1
						JMP		for4
						
						
				fine4:	SHL		EBX,2
						MOV	ESI, [EBP+12]			;riprendo la size
						
						
				fors4:	CMP		EBX,ESI
						JNL		return4
						SHL		EBX,2
						MOVSS	XMM1,[ECX+EBX]		;g[h]
						ADD		EBX,EDI
						;MOVUPS	XMM2,[EDX+EDI]		;G[riga*size+h]
						SQRTSS	XMM2,[EDX+EBX]		;radice G[riga*size+h]
						SUB		EBX,EDI
						ADDSS	XMM2,XMM0			;radice+epsilon
						DIVSS	XMM1,XMM2
						MOVSS	XMM3,[EAX+EBX]
						ADDSS	XMM3,XMM1
						MOVSS	[EAX+EBX],XMM3
						
						SHR		EBX,2
						ADD		EBX,1
						JMP		fors4
						
				
				return4:	
						pop	edi									; ripristina i registri da preservare
						pop	esi
						pop	ebx
						mov	esp, ebp								; ripristina lo Stack Pointer
						pop	ebp									; ripristina il Base Pointer
						ret										; torna alla funzione C chiamante

global calcoloTheta

calcoloTheta:
				
				push		ebp							; salva il Base Pointer
				mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
				push		ebx							; salva i registri da preservare
				push		esi
				push		edi
				
				
				MOV	EAX,[EBP+8]		;vector daSottrarre
				MOV	EBX,[EBP+12]		;vector Theta
				MOVSS	XMM0,[EBP+16]	;eta/v
				MOV	ECX,[EBP+20]		;size
				
				SHUFPS	XMM0,XMM0,00000000
				
				XOR		ESI,ESI
				MOV	EDI,ECX			;copia di size
		
				SHR		ECX,3			;size//4
				
				for5:		CMP		ESI,ECX
						JNL		fine5
						SHL		ESI,5				;h*dim
						
						MOVUPS	XMM1,[EAX+ESI]		;DaSottrarre[h]
						MOVUPS	XMM2,[EAX+ESI+16]	
						MOVUPS	XMM3,[EAX+ESI+32]	
						;MOVUPS	XMM4,[EAX+ESI+48]	
						
						MULPS	XMM1,XMM0			;daSottrarre[h]*eta/v
						MULPS	XMM2,XMM0
						MULPS	XMM3,XMM0
						;MULPS	XMM4,XMM0
						
						;MOVUPS	[EAX+ESI],XMM1		;non lo inserisco in memoria perchè verrà azzerato
						MOVUPS	XMM4,[EBX+ESI]		;theta[h]
						MOVUPS	XMM5,[EBX+ESI+16]
						MOVUPS	XMM6,[EBX+ESI+32]
						;MOVUPS	XMM0,[EBX+ESI+48]
						
						SUBPS	XMM4,XMM1
						SUBPS	XMM5,XMM2
						SUBPS	XMM6,XMM3
						;SUBPS	XMM0,XMM4
						
						MOVUPS	[EBX+ESI],XMM4
						MOVUPS	[EBX+ESI+16],XMM5
						MOVUPS	[EBX+ESI+32],XMM6
						;MOVUPS	[EBX+ESI+48],XMM0
						
						SHR		ESI,5
						ADD		ESI,1
						JMP		for5
						
						
				fine5:	SHL		ECX,3
						
						
				fors5:	CMP		ECX,EDI
						JNL		return5
						SHL		ECX,2				;h*dim
						MOVSS	XMM1,[EAX+ECX]		;DaSottrarre[h]
						MULSS	XMM1,XMM0			;daSottrarre[h]*eta/v
						;MOVUPS	[EAX+ECX],XMM1		;non lo inserisco in memoria perchè verrà azzerato
						MOVSS	XMM2,[EBX+ECX]		;theta[h]
						SUBSS	XMM2,XMM1
						MOVSS	[EBX+ECX],XMM2
						SHR		ECX,2
						ADD		ECX,1
						JMP		fors5
						
						
				
				return5:	
						pop	edi									; ripristina i registri da preservare
						pop	esi
						pop	ebx
						mov	esp, ebp								; ripristina lo Stack Pointer
						pop	ebp									; ripristina il Base Pointer
						ret										; torna alla funzione C chiamante

global azzeramentoDaSottrarre
				
azzeramentoDaSottrarre:

				push		ebp							; salva il Base Pointer
				mov			ebp, esp					; il Base Pointer punta al Record di Attivazione corrente
				push		ebx							; salva i registri da preservare
				push		esi
				push		edi
				
				MOV		EAX,[EBP+8]		;daSottrarre
				MOV		EBX,[EBP+12]		;size
				
				MOV		ECX,EBX			;copia di size
				XOR			ESI,ESI
				
				SHR			EBX,2			;size//4
				XORPS		XMM0,XMM0
				
				
				for6:	CMP	 	ESI, EBX
					JNL		 fine6
					SHL		ESI, 4				;i*16
					
					
					;NON CONVIENE FARE LOOP UNROLLING PERCHE' CI SONO SOLO ACCESSI IN MEMORIA
					MOVUPS	[EAX+ESI],XMM0
					;MOVUPS	[EAX+ESI+16],XMM0
					;MOVUPS	[EAX+ESI+32],XMM0
					;MOVUPS	[EAX+ESI+48],XMM0
					
					SHR		ESI, 4				;i/16
					ADD		ESI, 1				;i++
					JMP	 	for6
					
				fine6: SHR	EBX,2
				
				fors6: CMP	ESI,ECX
					  JNL	return6
					  SHL	ESI,2
					  
					  MOVSS	[EAX+ESI],XMM0
					  
					  SHR	ESI,2
					  ADD	ESI,1
					  JMP	fors6
					  
				return6:	pop	edi									; ripristina i registri da preservare
						pop	esi
						pop	ebx
						mov	esp, ebp							; ripristina lo Stack Pointer
						pop	ebp									; ripristina il Base Pointer
						ret										; torna alla funzione C chiamante
				
