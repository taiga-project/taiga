pro am_renate110taiga;,atom_,n_point,E_beam,endpoints_main,angle, atom_=atom_,E_beam=E_beam,$
;endpoints_main=endpoints_main,angle=angle

args = command_line_args(count=args_n)
;for i=0,args_n-1 do begin
;  print,args[i]
;endfor

default,n_point,200
default,angle,0
default,E_beam,60
default,atom_,'Li'
default,shotnumber,'30866';'11901'
default,appliance,'compass'
default,time_,'106';'1150'
default,deflection_plane_pos,2.30
default,endpoints_main,[[0.78,0.0],[0.35,0.0]]


;DEFAULT INPUT
if args_n gt 0 then begin 
  shotnumber=args[0]
  time_=args[1]
  atom_=args[2]
  E_beam=float(args[3])
  angle=float(args[4])  
  print,'Input OK.'
end

if args_n gt 5 then begin  
  n_point=float(args[5])
  print,'Npoint input OK.'
end

if args_n gt 6 then begin  
  deflection_plane_pos=float(args[6])
  print,'Deflection plane position input OK.'
end

if args_n gt 7 then begin  
  appliance=args[7]
  print,'Appliance input OK.'
end

if args_n gt 8 then begin 
  endpoints_main[0,0]=float(args[8])
  endpoints_main[1,0]=float(args[9])
  endpoints_main[0,1]=float(args[10])
  endpoints_main[1,1]=float(args[11])  
  print,'Endppoint input OK.'
end
 print,endpoints_main[1,0]
 endpoints_main[1,0] += sin(angle/180.0*!pi)*(deflection_plane_pos-endpoints_main[0,0])

 print,endpoints_main[1,0]
 print,endpoints_main[*,*]
 
cd,current=c
print,c

pi_virtual10, appliance=appliance, shotnumber=shotnumber, time_=time_, n_beam=0, $
E_beam=E_beam, atom_=atom_, endpoints_main, angle=angle, n_point=n_point,startlength=0,    /write_states;,/save_all
return
end

