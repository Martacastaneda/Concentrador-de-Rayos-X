import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
import matplotlib.pyplot as plt3

from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import random as pepe

# generate random floating point values
from random import seed
from random import random
from random import randint
from random import gauss
from random import choices

# seed random number generator
seed(1)
archivo_csv = '/Users/jruzarme/Desktop/CodigoMarta/nosee.csv'
nuevo_df = pd.read_csv(archivo_csv, delimiter=';')
print(nuevo_df)
nuevo_df.to_csv('datos.txt', sep='\t', index=False,  float_format='%.2f')

def obtener_reflectividad_desde_csv(archivo_csv, angulo, energia):
    energiab=energia[0]
    nuevo_df = pd.read_csv(archivo_csv, delimiter=';')
    idx = (nuevo_df.iloc[:, 0] - angulo).abs().idxmin()
    columns_without_unnamed = [col for col in nuevo_df.columns if 'Unnamed' not in col]
    closest_energy_col = min(columns_without_unnamed, key=lambda x: abs(float(x) - energiab))
    reflectividad = nuevo_df.at[idx, closest_energy_col]

    if not pd.isnull(reflectividad):
        return reflectividad
    else:
        return 0.0  # Valor predeterminado si no hay datos de reflectividad


def dibujar_rectangulos():
    # Crear una figura y ejes para el ray-tracing
    fig,ax  = plt.subplots()
    # Configurar los límites del gráfico
    ax.set_xlim(-100, 100)  # Ajuste los límites según tus necesidades
    ax.set_ylim(-20, 20)  # Ajuste los límites según tus necesidades
        
    
    # Especificar las coordenadas de la fuente 
    fuente_x = 0  # Coordenada x del centro del primer rectángulo
    fuente_y = 0  # Coordenada y del centro del primer rectángulo
    fuente_width = 20  # Ancho del primer rectángulo 
    fuente_height = 55  # Altura total del primer rectángulo

    # Crear el primer rectángulo
    fuente = Rectangle((fuente_x - fuente_width/2, fuente_y - fuente_height/2), fuente_width, fuente_height, facecolor='blue')

    # Agregar el primer rectángulo al gráfico
    ax.add_patch(fuente)


    # Especificar las coordenadas del detector
    detector_x = fuente_x + 601 #3360+3 #1781  # 1000  # Coordenada x del centro del tercer rectángulo
    detector_y = 0  # Coordenada y del centro del tercer rectángulo
    detector_width = 20  # Ancho del tercer rectángulo (ajusta según tus necesidades)
    detector_height = 150  # Altura del tercer rectángulo

    # Crear el tercer rectángulo
    detector = Rectangle((detector_x - detector_width/2, detector_y - detector_height/2), detector_width, detector_height, facecolor='black')

    # Agregar el tercer rectángulo al gráfico
    ax.add_patch(detector)




    # Especificar las coordenadas del espejo frontal-superior  <<<-----!!!
    rect2_x = fuente_x + 300  # Coordenada x del centro del segundo rectángulo
    rect2_y = 15  # Coordenada y del centro del segundo rectángulo
    rect2_width = 300  # Ancho del segundo rectángulo 
    rect2_height = 5  # Altura del segundo rectángulo
    angle_A=2.0## aqui el angulo esrta bien para el rectángulo

    # Crear el segundo rectángulo
    rect2 = Rectangle((rect2_x - rect2_width/2, fuente_y+fuente_height/2), rect2_width, rect2_height, facecolor='red', angle=-angle_A)

    # Agregar el segundo rectángulo al gráfico
    ax.add_patch(rect2)

    # Especificar las coordenadas del espejo frontal-inferior   <<<-----!!!
    rect3_x = fuente_x + 300  # Coordenada x del centro del segundo rectángulo
    rect3_height = 5  # Altura del segundo rectángulo
    rect3_y = rect2_y# Coordenada y del centro del segundo rectángulo
    rect3_width = 300  # Ancho del segundo rectángulo
    
    #angle_A=15 ## aqui el angulo esrta bien para el rectángulo

    # Crear el segundo rectángulo
    rect3 = Rectangle((rect3_x - rect3_width/2, fuente_y-fuente_height/2), rect3_width, -rect3_height, facecolor='red', angle=angle_A)

    # Agregar el segundo rectángulo al gráfico
    ax.add_patch(rect3)



    # Cálculo ecuación del plano del espejo frontal-superior
    #
    x1=rect2_x-rect2_width/2
    y1=fuente_y+fuente_height/2
    x2=x1+rect2_width*np.cos(np.radians(angle_A))
    y2=y1-rect2_width*np.sin(np.radians(angle_A))
    p_term=(y1-y2)/(x1-x2)
    b_term=y1-p_term*x1
    print()
    print("E1 ::", "Pendiente",p_term, "Independiente",b_term)

    plt.plot([x1,x2],[y1,y2],color='blue', linestyle=':', linewidth=1)
      
      
    # Cálculo ecuación del plano del espejo frotal-inferior
    #
    x3=rect3_x-rect3_width/2-rect3_height*np.sin(np.radians(angle_A))
    y3=fuente_y-fuente_height/2+ rect3_height*np.cos(np.radians(angle_A))-rect3_height
    x4=x3+rect3_width*np.cos(np.radians(-angle_A))
    y4=y3+rect3_width*np.sin(np.radians(angle_A))
    p2_term=(y3-y4)/(x3-x4)
    b2_term=y3-p2_term*x3
    print()
    print("E2 ::", "Pendiente",p2_term, "Independiente",b2_term)

    plt.plot([x3,x4],[y3,y4],color='blue', linestyle=':', linewidth=1)
    
    ########### Hasta aquí solo dibujamos rectángulos y las líneas de corte de los espejos


    # Crear una lista para almacenar los datos de la tabla
    tabla_datos = []

    #############################################################(JRA_1)
    #############################################################
    # Comienzon sampleador de espectro de axiones
    #
    # Axion flux [cm−2 s−1 keV−1]:  6.02 × 10^10  g10^2 E^{2.481}exp{−E/1.205}
    #
    def f(E):
        a=6.02e10
        return a*np.power(E,2.481)*np.exp(-E/1.205)
    
    x_range = np.linspace(0,12, num=24)
    integral = np.trapz(f(x_range),x_range)
   
    def w(E):
        b=(1/integral)*0.5
        return b*f(E)
    
    print("Wow integramos [0-to-12 keV]:", integral)
    prueba =(1/integral)*(0.5)*(f(0.5)+f(1)+f(1.5)+f(2)+f(2.5)+f(3)+f(3.5)+f(4)+f(4.5)+f(5)+f(5.5)+f(6)+f(6.5)+f(7)+f(7.5)+f(8)+f(8.5)+f(9)+f(9.5)+f(10)+f(10.5)+f(11)+f(11.5)+f(12)  )
    print("Test normalization", prueba)
    
    testList   = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
                  6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0]
     
    weights =[w(0.0),w(0.5), w(1.0), w(1.5), w(2.0), w(2.5), w(3.0), w(3.5), w(4.0), w(4.5), w(5.0), w(5.5), w(6.0),
                  w(6.5), w(7.0), w(7.5), w(8.0), w(8.5), w(9.0), w(9.5), w(10.0), w(10.5), w(11.0), w(11.5), w(12.0)]
                  
    #print("Energy photon [keV]:", testList)
    #print("Energy probability [1]:",weights)

    

    #  Final del sampleador del espectro de axiones pesado (obtienes mas puntos en funcion de la forma del espectro)
    #
    #  Para cada rayo tenemos que elegir su energía. Se hace con el comando choices:
    #
    #  randomSampleList = choices(testList, weights)
    ##############################
    ##############################(JRA_1f)

    rounds=1 #rounds=1 makes one revolution    x4 makes lots of statistics rounds
    #for rotate in range(0, 1, 1):
    for rotate in range(0, rounds*180, 1):
        #angle=rotate
        #angle=randint(0,180)
        angle=randint(0,180) + gauss(0, 1)
        # Actualizar las coordenadas de los espejos con el ángulo actual
        #angle_rad = np.radians(angle)
        rect2.set_xy((rect2_x - rect2_width/2, fuente_y + fuente_height/2))
        rect2.set_angle(-angle_A)
        rect3.set_xy((rect3_x - rect3_width/2, fuente_y - fuente_height/2))
        rect3.set_angle(angle_A)

        # Calcular la distancia entre los rayos (100 micron)
        distancia_rayos = 0.5  #was 0.1

        # Calcular el número de rayos a trazar
        num_samples = int(fuente_height / distancia_rayos)
        #num_samples =1100
        
        #Escribimos angulo simulado, distancia rayos y numero de rayos simulados
        #print("angle ::", angle, "dist_ray ::", distancia_rayos, "num_ray ::", num_samples)
        print("[",rotate,"]  ","angle ::", angle, "num_ray ::", num_samples)


        for i in range(num_samples):
            # Calcular la coordenada y del rayo actual
            #distancia_rayos = random()
            rayo_x = fuente_x + fuente_width / 2
            rayo_y = fuente_y + fuente_height/2 - i * distancia_rayos
            #rayo_y = randint(0,fuente_height)+distancia_rayos


            # Calcular la coordenada x final en la pantalla si choca con cualquier superficie del espejo
            punto_final_x = detector_x
            #print("Xinicial",rayo_y, "Xdetector",punto_final_x)


            #############################################(JRA_2)
            #############################################
            #### Esto selecciona la energía del fotón
            ##
            ##
            EnergiaFoton = choices(testList, weights)
            #print("Energia Foton [keV]: ",EnergiaFoton)
            #
            # Fin de selección de energía
            ##############################################(JRA2_f)
            

            ######################
            ###  Toca cambiar esto
            ##
            for jj in range(-135,135,1):
                
                #Calcular corte con espejo E1(frontal-superior) ... f(x)_E1 = rayo_y
                
                
                #Para emision colimada
                #xE1= (rayo_y-b_term)/p_term
                #yE1=rayo_y
                #Para fuente extensa plana
                if((np.tan(np.radians(jj/3))-p_term)!=0):
                    xE1= b_term/(np.tan(np.radians(jj/3))-p_term)  ###careful here with division
                    yE1=np.tan(np.radians(jj/3))*xE1
                else:
                    xE1=None
                    yE1=None
                    reflectividad=0.0
                #Calcular corte con espejo E2(frontal-inferior) ... f(x)_E1 = rayo_y
                
                #Para emision colimada
                #xE2= (rayo_y-b2_term)/p2_term
                #yE2=rayo_y
                #Para fuente extensa plana
                if ((np.tan(np.radians(jj/3))-p2_term)!=0):
                    xE2= b2_term/(np.tan(np.radians(jj/3))-p2_term)   ####careful here with division
                    yE2=np.tan(np.radians(jj/3))*xE2
                else:
                    xE2=None
                    yE2=None
                    reflectividad=0.0
                
                #if jj==-45 or jj==45:
                #    print("[",jj,"] test")

                ####Vale, esto está funcionando ....
                #Condicion de estar en el Espejo Frontal superior <------!!!!
                if (yE1!=None) and (yE1<=y1 and yE1>=y2):
                    plt.plot(xE1,yE1,marker="o", markersize=0.5, markeredgecolor="yellow")
                    slope_out = 2*p_term
                    new_b_term=yE1-slope_out*xE1
                    x_final=detector_x
                    y_final=slope_out*x_final+new_b_term

            
                else:
                    # Rayo no choca con el espejo
                    xE1 = None
                    reflectividad=0.0
        
                if xE1 is not None:
                    # Trazar el rayo original hasta el punto de impacto
                    outray_x = [xE1, x_final]
                    outray_y = [yE1, y_final]
                    ####ax.plot(outray_x, outray_y, color='limegreen', linestyle=':', linewidth=1)
    
                    inray_x = [rayo_x, xE1]
                    inray_y = [rayo_y, yE1]
                    ####ax.plot(inray_x, inray_y, color='magenta', linestyle=':', linewidth=1)
            
                    #ax2.plot2(x_final,y_final,marker="o", markersize=0.5, markeredgecolor="yellow")
                    # Calcular el ángulo con el que el rayo choca con el espejo (ángulo de incidencia)
                    pendiente_incidencia = (yE1 - rayo_y) / (xE1 - rayo_x)
                    ##--> angulo_choque = np.degrees(np.arctan(pendiente_incidencia))  ## esto estaba mal, hay que añadir angulo espejo
                    angulo_choque = np.degrees(np.arctan(pendiente_incidencia))-np.degrees(np.arctan(p_term))
                   
                    ##print("[",i,"::",jj,"]  ","graze_angle ::", angulo_choque)

                   
                    # Calcular el ángulo con el que el rayo sale del espejo (ángulo de reflexión)
                    pendiente_reflexion = (y_final - yE1) / (x_final - xE1)
                    angulo_salida = np.degrees(np.arctan(pendiente_reflexion))

                    # Calcular reflectividad usando la función modificada
                    #reflectividad = obtener_reflectividad_desde_csv(archivo_csv, abs(angulo_choque), EnergiaFoton[0])
                    reflectividad = obtener_reflectividad_desde_csv(archivo_csv, abs(angulo_choque), EnergiaFoton)

                    if abs(angulo_choque) > 15.5:
                        reflectividad = 0.0



                    # Agregar los datos a la tabla
                    datos = {
                        'Ángulo': angle,
                        'Posición Fuente (xF)': fuente_x+fuente_width/2,
                        'Posición Fuente (yF)': rayo_y,
                        'Posición Impacto (x_FrontSup)': xE1,
                        'Posición Impacto (y_FrontSup)': yE1,
                        'Posición Final (xD)': x_final,
                        'Posición Final (yD)': y_final,
                        'Ángulo Choque (grados)': angulo_choque,
                        'Ángulo Salida (grados)': angulo_salida,
                        'Reflectividad' : reflectividad,
                        'Energía Fotón': float(EnergiaFoton[0])
                        
                    }
                    tabla_datos.append(datos)
        
        
        
                else:
                    # Agregar los datos a la tabla
                    """
                    reflectividad=0
                    datos = {
                        'Ángulo': angle,
                        'Posición Fuente (xF)': fuente_x+fuente_width/2,
                        'Posición Fuente (yF)': rayo_y,
                        'Posición Impacto (x_FrontSup)': "-",
                        'Posición Impacto (y_FrontSup)': "-",
                        'Posición Final (xD)': "-",
                        'Posición Final (yD)': "-"
                        
                    }
                    tabla_datos.append(datos)"""


                ####Vale, esto está funcionando ....
                #Condicion de estar en el Espejo Frontal inferior <------!!!!
                if (yE2!=None) and (yE2>=y3 and yE2<=y4):
                    plt.plot(xE2,yE2,marker="o", markersize=0.5, markeredgecolor="yellow")
        
                    slope_out2 = 2*p2_term
                    new_b_term2=yE2-slope_out2*xE2
                    x_final2=detector_x
                    y_final2=slope_out2*x_final2+new_b_term2

            
                else:
                    # Rayo no choca con el espejo
                    xE2 = None
                    reflectividad=0.0


                #Writing data to table
                if xE2 is not None:
                    # Trazar el rayo original hasta el punto de impacto
                    outray_x2 = [xE2, x_final2]
                    outray_y2 = [yE2, y_final2]
                    ####ax.plot(outray_x2, outray_y2, color='limegreen', linestyle=':', linewidth=1)

                    inray_x2 = [fuente_width/2, xE2]
                    inray_y2 = [rayo_y, yE2]
                    ####ax.plot(inray_x2, inray_y2, color='magenta', linestyle=':', linewidth=1)
            
                    #ax2.plot2(x_final2,y_final2,marker="o", markersize=0.5, markeredgecolor="yellow")
                    # Calcular el ángulo con el que el rayo choca con el espejo E2 (ángulo de incidencia)
                    pendiente_incidencia_E2 = (yE2 - rayo_y) / (xE2 - (fuente_x + fuente_width / 2))
                    angulo_choque_E2 = -np.degrees(np.arctan(pendiente_incidencia_E2))+np.degrees(np.arctan(p2_term))

                    # Calcular el ángulo con el que el rayo sale del espejo E2 (ángulo de reflexión)
                    pendiente_reflexion_E2 = (y_final2 - yE2) / (x_final2 - xE2)
                    angulo_salida_E2 = np.degrees(np.arctan(pendiente_reflexion_E2))

                    #reflectividad_E2 = obtener_reflectividad_desde_csv(archivo_csv, abs(angulo_choque_E2), EnergiaFoton[0])
                    reflectividad_E2 = obtener_reflectividad_desde_csv(archivo_csv, abs(angulo_choque_E2), EnergiaFoton)

                    if abs(angulo_choque_E2) > 15.5:
                        reflectividad_E2 = 0.0


                    # Agregar los datos a la tabla
                    datos = {
                        'Ángulo': angle,
                        'Posición Fuente (xF)': fuente_x+fuente_width/2,
                        'Posición Fuente (yF)': rayo_y,
                        'Posición Impacto (x_FrontInf)': xE2,
                        'Posición Impacto (y_FrontInf)': yE2,
                        'Posición Final (xD)': x_final2,
                        'Posición Final (yD)': y_final2,
                        'Ángulo Choque (grados)': angulo_choque_E2,
                        'Ángulo Salida (grados)': angulo_salida_E2,
                        'Reflectividad' : reflectividad_E2,
                        'Energía Fotón': float(EnergiaFoton[0])
                    }
                    tabla_datos.append(datos)
        
        
                else:
                    # Agregar los datos a la tabla
                    """
                    datos = {
                        'Ángulo': angle,
                        'Posición Fuente (xF)': fuente_x+fuente_width/2,
                        'Posición Fuente (yF)': rayo_y,
                        'Posición Impacto (x_FrontInf)': "-",
                        'Posición Impacto (y_FrontInf)': "-",
                        'Posición Final (xD)': "-",
                        'Posición Final (yD)': "-"
                    }
                    tabla_datos.append(datos)"""




    # Mostrar el gráfico
    plt.axis('equal')  # Para mantener las proporciones de los ejes
    plt.show()

    # Crear un DataFrame a partir de los datos de la tabla
    df = pd.DataFrame(tabla_datos)
    df.to_csv('datos_rayos_fuente_noo_colimada.txt', sep='\t', index=False)
    print(df)
    return df 




def dibujar_hits(df):
    # Crear una figura y ejes
    fig2,ax2 = plt.subplots()
    
    # Configurar los límites del gráfico
    ax2.set_xlim(-50, +50)  # Ajuste los límites según tus necesidades
    ax2.set_ylim(-50, 50)  # Ajuste los límites según tus necesidades
    # Dibujar los puntos de llegada de los rayos al detector
    
    #ax2.scatter(df['Ángulo'], df['Posición Final (yD)'], color='red', marker=',', s=1)
    #Aquí usamos una matriz de rotación en función del ángulo de simulación:
    #(x1)    ( cos   sen)(x0)
    #    =
    #(y1)    ( -sen  cos)(y0)
    #
    #Cuidado, porque simulación tiene siempre x0=0
    #La rotación queda:
    #x1=sen .y0
    #y1=cos .y0
    #
    ax2.scatter(np.sin(np.radians(df['Ángulo']))*df['Posición Final (yD)'], np.cos(np.radians(df['Ángulo']))*df['Posición Final (yD)'], color='red', marker='.', s=1)
    
    # Agregar etiquetas con valores numéricos en los puntos
    for index, row in df.iterrows():
        x = row['Posición Final (xD)']
        y = row['Posición Final (yD)']
        #valor = f"({x:.2f}, {y:.2f})"  # Formatea el valor numérico como (x, y) con dos decimales
        #ax2.annotate(valor, (x, y), textcoords="offset points", xytext=(0, 10), ha='center')
    plt.show()
    """

    for index, row in df.iterrows():
        xE1 = row['Posición Impacto (x_FrontSup)']
        yE1 = row['Posición Impacto (y_FrontSup)']
        xE2 = row['Posición Impacto (x_FrontInf)']
        yE2 = row['Posición Impacto (y_FrontInf)']
        x_final = row['Posición Final (xD)']
        y_final = row['Posición Final (yD)']

        # Verifica si el punto de corte con el espejo frontal superior es válido (no es None)
        if xE1 is not None:
            ax2.scatter(xE1, yE1, color='green', marker='o', s=20)
            ax2.plot([xE1, x_final], [yE1, y_final], color='green', linestyle=':', linewidth=1)

        # Verifica si el punto de corte con el espejo frontal inferior es válido (no es None)
        if xE2 is not None:
            ax2.scatter(xE2, yE2, color='blue', marker='o', s=20)
            ax2.plot([xE2, x_final], [yE2, y_final], color='blue', linestyle=':', linewidth=1)

        # Dibujar puntos de llegada al detector
        ax2.scatter(x_final, y_final, color='red', marker='o', s=50)
    plt.show"""
    
      
    """ 
    # Especificar las coordenadas de la fuente
    bla_x = 0  # Coordenada x del centro del primer rectángulo
    bla_y = 0  # Coordenada y del centro del primer rectángulo
    bla_width = 150  # Ancho del primer rectángulo
    bla_height = 150  # Altura total del primer rectángulo

    # Crear el primer rectángulo
    detector_plane = Rectangle((bla_x - bla_width/2, bla_y - bla_height/2), bla_width, bla_height, facecolor='blue')  
    # Agregar el primer rectángulo al gráfico
    ax2.add_patch(detector_plane)"""

    




# Ejecutar la función para dibujar los rectángulos y trazar los rayos de luz
df=dibujar_rectangulos()
dibujar_hits(df)

# Calcular el histograma de una columna específica, por ejemplo, 'Posición Final (xD)'
hist_data = df['Posición Final (yD)']

# Crear una nueva figura y ejes para el histograma
fig3, ax3 = plt.subplots()

# Trazar el histograma con la función plt.hist()
ax3.hist(hist_data, bins=100, color='skyblue', edgecolor='black')  # Puedes ajustar el número de bins según tus necesidades

# Agregar etiquetas y título al histograma
ax3.set_xlabel('Valores')
ax3.set_ylabel('Frecuencia')
ax3.set_title('Histograma de Valores')
ax3.set_xticks(range(-40, 40, 5))
#ax3.set_yticks(range(0, 9, 1))

# Mostrar la tercera figura
plt.show()

###################################################(JRA_3)
###################################################
# Hacemos un plot al final para comprobar
# Aquí podíamos usar df['Energía Fotón'], sería más realista, pero eso lo haremos con el run de gran estadística
#
#

fig,ax_new  = plt.subplots()

def f(E):
    a=6.02e10
    return a*np.power(E,2.481)*np.exp(-E/1.205)
    
x_range = np.linspace(0,12, num=24)
integral = np.trapz(f(x_range),x_range)
   
def w(E):
    b=(1/integral)*0.5
    return b*f(E)
    

testList2   = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, 11.5, 12.0]
     
weights2 =[w(0.0),w(0.5), w(1.0), w(1.5), w(2.0), w(2.5), w(3.0), w(3.5), w(4.0), w(4.5), w(5.0), w(5.5), w(6.0), w(6.5), w(7.0), w(7.5), w(8.0), w(8.5), w(9.0), w(9.5), w(10.0), w(10.5), w(11.0), w(11.5), w(12.0)]
            
            
rounds=1
angles=1*180
divergence=3*90
start_points=110
value=rounds*angles*divergence*start_points
                        
randomSampleList2 = choices(testList2, weights2,k=value)

testing=float(df.shape[0])
factor=testing/value
print("Rows ::",df.shape[0],"-->",testing,"[",factor,"]")

plt3.hist(randomSampleList2,ec="black", bins=24)
plt3.hist(randomSampleList2,weights=factor*np.ones_like(randomSampleList2),color="orange",ec="black", bins=24)
plt3.hist(df['Energía Fotón'],weights=df['Reflectividad'],color="mediumorchid",ec="black",bins=24)

ax_new.set_xlabel('Energy [keV]')
ax_new.set_ylabel('Frecuencia')
ax_new.set_title('Simulated Axion Spectra')
#ax_new.set_xticks(range(0, 12, 1))
    

plt3.show()

#
# Fin plot espectro axiones
############################################(JRA_3f)


import matplotlib.pyplot as plt

def histo2D_ponderado(df):

    # Obtener las coordenadas x e y del scatter plot desde el DataFrame ... Cuidado=, hay que tener en cuanta el ángulo
    x_hist2D      = np.sin(np.radians(df['Ángulo']))*df['Posición Final (yD)'] #df['Posición Final (xD)']
    y_hist2D      = np.cos(np.radians(df['Ángulo']))*df['Posición Final (yD)'] #df['Posición Final (yD)']
    reflectividad = 100*df['Reflectividad']

        # Crear una nueva figura y ejes para el scatter plot
    fig_hist2D, ax_hexbin = plt.subplots()

    ax_hexbin.set_xlim(-10, 10)  # Ajuste los límites según tus necesidades
    ax_hexbin.set_ylim(-10, 10)  # Ajuste los límites según tus necesidades
    

   
    # Trazar el scatter plot ponderado
    #scatter = ax_scatter.scatter(x_scatter, y_scatter,c=reflectividad, s=sizes, alpha=0.2, cmap='viridis')
    #histo2D = ax_hexbin.hexbin(x_hist2D, y_hist2D, gridsize=40, bins='log', cmap='Blues')
    histo2D = ax_hexbin.hexbin(x_hist2D, y_hist2D, C=reflectividad, gridsize=30, cmap='Blues')

    # Agregar una barra de color para visualizar los tamaños
    cbar = fig.colorbar(histo2D,label='hits per bin')
    #cbar.set_label('Hits')

    # Agregar etiquetas y título al scatter plot
    ax_hexbin.set_xlabel('Posición Final (xD)')
    ax_hexbin.set_ylabel('Posición Final (yD)')
    ax_hexbin.set_title('Scatter Plot Ponderado por Reflectividad')


    
    sums=df['Reflectividad'].sum()
    
    print("Total ::", sums)
    print("")


    # Mostrar el scatter plot
    plt.show()

# Llamar a la función para crear el scatter plot ponderado
histo2D_ponderado(df)


###################################################(JRA_4)
###################################################
# Hacemos un plot al final para comprobar el espectro reflejado
#
#
#
import matplotlib.pyplot as pltNew

def histo_RefAXION(df):

    figRef, ax_newAxion  = pltNew.subplots()


    print(df)

    #df['EnergyPhoton_Float'] = df['Energía Fotón'].astype(float)
    
    print(df.dtypes)

    
    pltNew.hist(df['Energía Fotón'],weights=df['Reflectividad'],color="mediumorchid",ec="black",bins=24)
    ax_newAxion.set_xlabel('Energy [keV]')
    ax_newAxion.set_ylabel('Frecuencia')
    ax_newAxion.set_title('Reflected Axion Spectra')
    ax_new.set_xticks(range(0, 12, 1))

    pltNew.show()
    
histo_RefAXION(df)
#
# Fin plot espectro axiones
############################################(JRA_4f)
