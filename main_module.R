setwd("C:/Users/zhyxi/Documents/R/digibog2d")

oxic_decay_base <- 0.035
anoxic_decay_base <- 0.00001
base_temp <- 6.29
Q10_oxic <- 3
Q10_anoxic <- 3
density <- 0.1
porosity <- 0.3
k_param_a <- 1586.8724
k_param_b <- 8
t_extent <- 55
annual_tsteps <- 525600
output_interval <- 1
x_extent <- 3
y_extent <- 27
spatial_step <- 200
pond_depth <- 0.25
mean_temp <- 11.35
t_step_multi <- 0.5
max_tstep <- 5

z_extent <- t_extent+2

col_wt_depth <- col_wt_depth_av <- col_wt_depth_av_su <- col_wt_sum <- col_wt_sum_su <- col_mass_per_area <-
  base_altitude <- water_table <- activation_status <- array(0,dim=c(x_extent,y_extent))
water_change <- wk_mean <- array(1,dim=c(x_extent,y_extent))
wet_proportion <- array(0,dim=c(x_extent,y_extent,z_extent))
layer_mass <- layer_attributes <- array(0,dim=c(x_extent,y_extent,z_extent,3))
transmissivity <- layer_storage <- array(1,dim=c(x_extent,y_extent,z_extent,2))

year_counter <- 1
week_counter <- 0
annual_tsteps <- 1
week_tsteps <- annual_tsteps/52
timestep <- week_tsteps
t_step_sum <- 0
t_step_sum_su <- 0
max_tstep <- max_tstep/525600
mean_t_step <- 0
sub_year_counter <- 0
sub_year_counter_su <- 0
output_counter <- 0
oxic_decay <- 0
anoxic_decay <- 0
x <- 1
y <- 1
z <- 1
no_layers <- 3
new_week <- 1

column_status <- as.vector(t(read.table("040_DigiBog_BB_IN_column_status.txt")))
baltitude <- as.vector(t(read.table("050_DigiBog_BB_IN_baltitude.txt")))
grab <- 1
for (x in 1:x_extent){
  for (y in 1:y_extent){
    activation_status[x,y] <- column_status[grab]
    base_altitude[x,y] <- baltitude[grab]
    grab <- grab+1
  }
}

for (x in 1:x_extent){
  for (y in 1:y_extent){
    if (activation_status[x,y]=="on"){
      layer_attributes[x,y,no_layers,1] <- pond_depth
      layer_attributes[x,y,no_layers,2] <- 0
      layer_attributes[x,y,no_layers,3] <- 1
    }
  }
}

for (x in 1:x_extent){
  for (y in 1:y_extent){
    if (activation_status[x,y]=="on"){
      layer_attributes[x,y,1,1] <- 2
      layer_attributes[x,y,1,2] <- 3153.6
      layer_attributes[x,y,1,3] <- 0.3
      transmissivity[x,y,1,1] <- layer_attributes[x,y,1,1]
      transmissivity[x,y,1,2] <- transmissivity[x,y,1,1]*layer_attributes[x,y,1,2]
      wet_proportion[x,y,1] <- 1
      water_table[x,y] <- transmissivity[x,y,1,1]
    }
  }
}

for (x in 1:x_extent){
  for (y in 1:y_extent){
    if (activation_status[x,y]=="on"){
      layer_mass[x,y,2,1] <- (9.3^2)*0.0001
      layer_mass[x,y,2,2] <- layer_mass[x,y,2,1]
      layer_mass[x,y,2,3] <- 1
      layer_attributes[x,y,2,1] <- layer_mass[x,y,2,1]/density
      layer_attributes[x,y,2,2] <- k_param_a*exp(k_param_b)
      layer_attributes[x,y,2,3] <- porosity
      transmissivity[x,y,2,1] <- transmissivity[x,y,1,1]+layer_attributes[x,y,2,1]
      transmissivity[x,y,2,2] <- transmissivity[x,y,1,1]+layer_attributes[x,y,2,1]*layer_attributes[x,y,2,2]
      wet_proportion[x,y,2] <- 1
      col_mass_per_area[x,y] <- layer_mass[x,y,2,1]
      water_table[x,y] <- transmissivity[x,y,2,1]
    }
  }
}

for (x in 2:(x_extent-1)){
  for (y in 1:(y_extent-1)){
    if (activation_status[x,y]=="diri"){
      water_table[x,y] <- transmissivity[x,y+1,2,1]
    }
  }
}

net_rain <- as.vector(t(read.table("020_DigiBog_BB_IN_net_rain.txt")))
temp <- as.vector(t(read.table("030_DigiBog_BB_IN_temp.txt")))
grab <- 1
t_step_sum_output <- col_t_step_output <- column_height_output <- wt_height_output <- mass_area_output <- 
  wt_depth_output <- wt_depth_summer_output <- layer_mass_output <- transmiss_output <- layer_wet_prop_output <- 
  k_profile_output <- vector()
repeat{
  print(paste("Model year",year_counter))
  repeat{
    if (new_week==1){
      rainfall <- net_rain[grab]
      temperature <- temp[grab]
      grab <- grab+1
      oxic_decay <- oxic_decay_base*Q10_oxic^((temperature-base_temp)/10)
      anoxic_decay <- anoxic_decay_base*Q10_anoxic^((temperature-base_temp)/10)
      t_step_sum_output <- c(t_step_sum_output,t_step_sum)
      t_step_sum <- 0
      t_step_sum_su <- 0
      new_week <- 0
    }
    
    sub_year_counter <- sub_year_counter+1
    if ((week_counter>=17)&(week_counter<=31)){
      sub_year_counter_su <- sub_year_counter_su+1
    }

    col_mass_per_area <- array(0,dim=c(x_extent,y_extent))

    timestep <- week_tsteps
    transmax <- 0

    for (x in 1:x_extent){
      for (y in 1:y_extent){
        if (activation_status[x,y]=="on"){
          layer_storage[x,y,1,1] <- layer_attributes[x,y,1,1]
          layer_storage[x,y,1,2] <- layer_attributes[x,y,1,1]*layer_attributes[x,y,1,3]
          for (z in 2:no_layers){
            layer_storage[x,y,z,1] <- layer_storage[x,y,z-1,1]+layer_attributes[x,y,z,1]
            layer_storage[x,y,z,2] <- layer_attributes[x,y,z,1]*layer_attributes[x,y,z,3]
          }
        }
      }
    }
    
    for (x in 1:x_extent){
      for (y in 1:y_extent){
        if (activation_status[x,y]=="on"){
          if (water_table[x,y] >= transmissivity[x,y,no_layers-1,1]){
            wk_mean[x,y] <- transmissivity[x,y,no_layers-1,2]
            if (wk_mean[x,y] > transmax) {transmax <- wk_mean[x,y]}
            wk_mean[x,y] <- wk_mean[x,y]/transmissivity[x,y,no_layers-1,1]
          } else if (water_table[x,y] < transmissivity[x,y,1,1]){
            wk_mean[x,y] <- layer_attributes[x,y,1,2]
            if (transmissivity[x,y,1,2] > transmax) {transmax <- transmissivity[x,y,1,2]}
            # wk_mean[x,y] <- wk_mean[x,y]/water_table[x,y]
          } else {
            for (z in 1:(no_layers-2)){
              if (water_table[x,y] == transmissivity[x,y,z,1]){
                wk_mean[x,y] <- transmissivity[x,y,z,2]
                if (wk_mean[x,y] > transmax) {transmax <- wk_mean[x,y]}
                wk_mean[x,y] <- wk_mean[x,y]/transmissivity[x,y,z,1]
              } else if (water_table[x,y] > transmissivity[x,y,z,1]){
                if (water_table[x,y] < transmissivity[x,y,z+1,1]){
                  wk_mean[x,y] <- transmissivity[x,y,z,2]+layer_attributes[x,y,z+1,2]*(water_table[x,y]-transmissivity[x,y,z,1])
                  if (wk_mean[x,y] > transmax) {transmax <- wk_mean[x,y]}
                  wk_mean[x,y] <- wk_mean[x,y]/water_table[x,y]
                }
              }
            }
          }
        }
      }
    }
    
    timestep <- 0.5*(0.5*spatial_step)^2*porosity/transmax*0.8
    timestep <- timestep * t_step_multi
    
    if (year_counter > 50){
      if (timestep > 50/525600){
        timestep <- 50/525600
      }
    } else if (year_counter <= 50){
      if (timestep > max_tstep){
        timestep <- max_tstep
      }
    }
    if ((timestep+t_step_sum) > week_tsteps){
      timestep <- week_tsteps-t_step_sum
      new_week <- 1
      week_counter <- week_counter+1
    }
    t_step_sum <- t_step_sum+timestep
    if ((week_counter>=17)&(week_counter<=31)){
      t_step_sum_su <- t_step_sum_su+timestep
    }
    
    x_flux <- y_flux <- array(0,dim=c(x_extent,y_extent))
    for (x in 1:(x_extent-1)){
      for (y in 2:(y_extent-1)){
        if ((activation_status[x,y]=="neu")&(activation_status[x+1,y]=="on")){
          x_flux[x,y] <- 0
          } else if ((activation_status[x,y]=="diri")&(activation_status[x+1,y]=="on")){
            x_flux[x,y] <- wk_mean[x+1,y]*((water_table[x,y]+water_table[x+1,y])/2)*
              (base_altitude[x,y]+water_table[x,y]-base_altitude[x+1,y]-water_table[x+1,y])*2*timestep
          } else if ((activation_status[x,y]=="on")&(activation_status[x+1,y]=="neu")){
            x_flux[x,y] <- 0
          } else if ((activation_status[x,y]=="on")&(activation_status[x+1,y]=="diri")){
            x_flux[x,y] <- wk_mean[x,y]*((water_table[x,y]+water_table[x+1,y])/2)*
              (base_altitude[x,y]+water_table[x,y]-base_altitude[x+1,y]-water_table[x+1,y])*2*timestep
          } else if ((activation_status[x,y]=="on")&(activation_status[x+1,y]=="on")){
          x_flux[x,y] <- (2*wk_mean[x,y]*wk_mean[x+1,y])/(wk_mean[x,y]+wk_mean[x+1,y])*((water_table[x,y]+water_table[x+1,y])/2)*
            (base_altitude[x,y]+water_table[x,y]-base_altitude[x+1,y]-water_table[x+1,y])*timestep
          }
      }
    }
    for (y in 1:(y_extent-1)){
      for (x in 2:(x_extent-1)){
        if ((activation_status[x,y]=="neu")&(activation_status[x,y+1]=="on")){
          y_flux[x,y] <- 0
          } else if ((activation_status[x,y]=="diri")&(activation_status[x,y+1]=="on")){
            y_flux[x,y] <- wk_mean[x,y+1]*((water_table[x,y]+water_table[x,y+1])/2)*
              (base_altitude[x,y]+water_table[x,y]-base_altitude[x,y+1]-water_table[x,y+1])*2*timestep
          } else if ((activation_status[x,y]=="on")&(activation_status[x,y+1]=="neu")){
            y_flux[x,y] <- 0
          } else if ((activation_status[x,y]=="on")&(activation_status[x,y+1]=="diri")){
            y_flux[x,y] <- wk_mean[x,y]*((water_table[x,y]+water_table[x,y+1])/2)*
              (base_altitude[x,y]+water_table[x,y]-base_altitude[x,y+1]-water_table[x,y+1])*2*timestep
          } else if ((activation_status[x,y]=="on")&(activation_status[x,y+1]=="on")){
            y_flux[x,y] <- (2*wk_mean[x,y]*wk_mean[x,y+1])/(wk_mean[x,y]+wk_mean[x,y+1])*((water_table[x,y]+water_table[x,y+1])/2)*
              (base_altitude[x,y]+water_table[x,y]-base_altitude[x,y+1]-water_table[x,y+1])*timestep
        }
      }
    }
    for (x in 2:(x_extent-1)){
      for (y in 2:(y_extent-1)){
        if (activation_status[x,y]=="on"){
          water_change[x,y] <- (x_flux[x-1,y]-x_flux[x,y]+y_flux[x,y-1]-y_flux[x,y])/(spatial_step^2)
          water_change[x,y] <- water_change[x,y]+rainfall*timestep
        }
      }
    }
    
    for (x in 2:(x_extent-1)){
      for (y in 2:(y_extent-1)){
        if (activation_status[x,y]=="on"){
          if (water_table[x,y] <= layer_storage[x,y,1,1]){
            marker <- 1
          } else {
            for (z in 1:(no_layers-1)){
              if (water_table[x,y] > layer_storage[x,y,z,1]){
                if (water_table[x,y] <= layer_storage[x,y,z+1,1]){
                  marker <- z+1
                } else if (water_table[x,y] > layer_storage[x,y,no_layers,1]){
                  marker <- no_layers
                  break
                }
              }
            }
          }
          if (water_change[x,y]>0){
            if (water_table[x,y] < layer_storage[x,y,no_layers,1]){
              z <- marker
              if (water_table[x,y]==layer_storage[x,y,z,1]){
                z <- z+1
              }
              repeat{
                test_depth <- layer_storage[x,y,z,2]*(layer_storage[x,y,z,1]-water_table[x,y])/layer_attributes[x,y,z,1]
                if (water_change[x,y] <= test_depth){
                  water_table[x,y] <- water_table[x,y]+water_change[x,y]/layer_attributes[x,y,z,3]
                  break
                } else {
                  if (z==no_layers){
                    water_table[x,y] <- layer_storage[x,y,z,1]
                    break
                  } else {
                    water_table[x,y] <- layer_storage[x,y,z,1]
                    water_change[x,y] <- water_change[x,y]-test_depth
                    z <- z+1
                  }
                }
              }
            }
          }
          if (water_change[x,y]<0){
            if (water_table[x,y]>0){
              z <- marker
              repeat{
                if (z==1){
                  test_depth <- layer_storage[x,y,1,2]*water_table[x,y]/layer_attributes[x,y,1,1]
                } else {
                  test_depth <- layer_storage[x,y,z,2]*(water_table[x,y]-layer_storage[x,y,z-1,1])/layer_attributes[x,y,z,1]
                }
                if (abs(water_change[x,y]) <= test_depth){
                  water_table[x,y] <- water_table[x,y]+water_change[x,y]/layer_attributes[x,y,z,3]
                  break
                } else {
                  if (z==1){
                    water_table[x,y] <- 0
                    break
                  } else {
                    water_table[x,y] <- layer_storage[x,y,z-1,1]
                    water_change[x,y] <- water_change[x,y]+test_depth
                    z <- z-1
                  }
                }
              }
            }
          }
        }
      }
    }
    
    for (x in 1:x_extent){
      for (y in 1:y_extent){
        if (activation_status[x,y]=="on"){
          col_wt_depth[x,y] <- transmissivity[x,y,no_layers-1,1]-water_table[x,y]
          col_wt_sum[x,y] <- col_wt_sum[x,y]+col_wt_depth[x,y]
          if ((week_counter >= 17)&(week_counter <= 31)){
            col_wt_sum_su[x,y] <- col_wt_sum_su[x,y]+col_wt_depth[x,y]
          }
        }
      }
    }
    
    for (x in 1:x_extent){
      for (y in 1:y_extent){
        if (activation_status[x,y]=="on"){
          for (z in 2:(no_layers-1)){
            layer_mass[x,y,z,1] <- layer_mass[x,y,z,1]*wet_proportion[x,y,z]*exp(-anoxic_decay*timestep*layer_mass[x,y,z,3])+
              layer_mass[x,y,z,1]*(1-wet_proportion[x,y,z])*exp(-oxic_decay*timestep)
            layer_mass[x,y,z,3] <- layer_mass[x,y,z,1]/layer_mass[x,y,z,2]
            layer_attributes[x,y,z,1] <- layer_mass[x,y,z,1]/density
            col_mass_per_area[x,y] <- col_mass_per_area[x,y]+layer_mass[x,y,z,1]
          }
        }
      }
    }
    
    for (x in 1:x_extent){
      for (y in 1:y_extent){
        if (activation_status[x,y]=="on"){
          for (z in 2:(no_layers-1)){
            layer_attributes[x,y,z,2] <- k_param_a*exp(k_param_b*layer_mass[x,y,z,3])
          }
        }
      }
    }
    
    for (x in 1:x_extent){
      for (y in 1:y_extent){
        if (activation_status[x,y]=="on"){
          transmissivity[x,y,1,1] <- layer_attributes[x,y,1,1]
          transmissivity[x,y,1,2] <- layer_attributes[x,y,1,1]*layer_attributes[x,y,1,2]
          for (z in 2:(no_layers-1)){
            transmissivity[x,y,z,1] <- transmissivity[x,y,z-1,1]+layer_attributes[x,y,z,1]
            transmissivity[x,y,z,2] <- transmissivity[x,y,z-1,2]+layer_attributes[x,y,z,1]*layer_attributes[x,y,z,2]
          }
        }
      }
    }
    
    for (x in 1:x_extent){
      for (y in 1:y_extent){
        if (activation_status[x,y]=="on"){
          for (z in 1:(no_layers-1)){
            if (z==1){
              if (water_table[x,y] >= transmissivity[x,y,z,1]){
                wet_proportion[x,y,1] <- 1
              } else {
                wet_proportion[x,y,1] <- water_table[x,y]/layer_attributes[x,y,1,1]
              }
            } else if (water_table[x,y] >= transmissivity[x,y,z,1]){
              wet_proportion[x,y,z] <- 1
            } else if (water_table[x,y] <= transmissivity[x,y,z-1,1]){
              wet_proportion[x,y,z] <- 0
            } else {
              wet_proportion[x,y,z] <- (water_table[x,y]-transmissivity[x,y,z-1,1])/layer_attributes[x,y,z,1]
            }
          }
        }
      }
    }
    
    if (week_counter==52) {break}
  }
  
  for (x in 1:x_extent){
    for (y in 1:y_extent){
      if (activation_status[x,y]=="on"){
        col_wt_depth_av[x,y] <- col_wt_sum[x,y]/sub_year_counter
        col_wt_depth_av_su[x,y] <- col_wt_sum_su[x,y]/sub_year_counter_su
      }
    }
  }
  
  mean_t_step <- 1/sub_year_counter*525600
  col_t_step_output <- c(col_t_step_output,mean_t_step)
  mean_t_step <- 0

  week_counter <- 0
  
  output_counter <- output_counter+1
  if (output_counter==output_interval){
    for (x in 1:x_extent){
      for (y in 1:y_extent){
        if (activation_status[x,y]=="on"){
          column_height_output <- c(column_height_output,transmissivity[x,y,no_layers-1,1])
          wt_height_output <- c(wt_height_output,water_table[x,y])
          mass_area_output <- c(mass_area_output,col_mass_per_area[x,y])
          wt_depth_output <- c(wt_depth_output,col_wt_depth_av[x,y])
          wt_depth_summer_output <- c(wt_depth_summer_output,col_wt_depth_av_su[x,y])
        } else {
          column_height_output <- c(column_height_output,-999)
          wt_height_output <- c(wt_height_output,-999)
          mass_area_output <- c(mass_area_output,-999)
          wt_depth_output <- c(wt_depth_output,-999)
          wt_depth_summer_output <- c(wt_depth_summer_output,-999)
        }
      }
    }
  }
  
  year_counter <- year_counter + 1
  
  if (year_counter > t_extent) {break}
  
  col_wt_sum <- array(0,dim=c(x_extent,y_extent))
  col_wt_sum_su <- array(0,dim=c(x_extent,y_extent))
  sub_year_counter <- 0
  sub_year_counter_su <- 0
  
  no_layers <- no_layers+1
  
  for (x in 1:x_extent){
    for (y in 1:y_extent){
      if (activation_status[x,y]=="on"){
        layer_attributes[x,y,no_layers,1] <- pond_depth
        layer_attributes[x,y,no_layers,2] <- 0
        layer_attributes[x,y,no_layers,3] <- 1
      }
    }
  }
  
  for (x in 1:x_extent){
    for (y in 1:y_extent){
      if (activation_status[x,y]=="on"){
        if (col_wt_depth_av[x,y] > 66.8){
          layer_mass[x,y,no_layers-1,1] <- 0.0000001
        } else {
          layer_mass[x,y,no_layers-1,1] <- (9.3+1.33*col_wt_depth_av[x,y]-
                                                 0.022*col_wt_depth_av[x,y]^2)^2*0.0001*(0.1575*mean_temp+0.009132)
        }
        layer_mass[x,y,no_layers-1,2] <- layer_mass[x,y,no_layers-1,1]
        layer_mass[x,y,no_layers-1,3] <- 1
        layer_attributes[x,y,no_layers-1,1] <- layer_mass[x,y,no_layers-1,1]/density
        layer_attributes[x,y,no_layers-1,3] <- porosity
        transmissivity[x,y,no_layers-1,1] <- transmissivity[x,y,no_layers-1-1,1]+layer_attributes[x,y,no_layers-1,1]
        if (water_table[x,y] >= transmissivity[x,y,no_layers-1,1]){
          wet_proportion[x,y,no_layers-1] <- 1
        } else if (water_table[x,y] <= transmissivity[x,y,no_layers-1,1]){
          wet_proportion[x,y,z] <- 0
        } else {
          wet_proportion[x,y,no_layers-1] <- (water_table[x,y]-transmissivity[x,y,no_layers-1,1])/layer_attributes[x,y,no_layers-1,1]
        }
        if (wet_proportion[x,y,no_layers-1] > 0){
          layer_attributes[x,y,no_layers-1,2] <- k_param_a*exp(k_param_b)
          transmissivity[x,y,no_layers-1,2] <- transmissivity[x,y,no_layers-1-1,2]+layer_attributes[x,y,no_layers-1,1]*layer_attributes[x,y,no_layers-1,2]
        }
        col_mass_per_area[x,y] <- col_mass_per_area[x,y]+layer_mass[x,y,no_layers-1,1]
      }
    }
  }
  
  output_counter <- 0
}

for (x in 1:x_extent){
  for (y in 1:y_extent){
    if (activation_status[x,y]=="on"){
      for (z in 2:(no_layers-1)){
        layer_mass_output <- c(layer_mass_output,layer_mass[x,y,z,3])
        transmiss_output <- c(transmiss_output,transmissivity[x,y,z,2])
        layer_wet_prop_output <- c(layer_wet_prop_output,wet_proportion[x,y,z])
        k_profile_output <- c(k_profile_output,layer_attributes[x,y,z,2])
      }
    }
  }
}