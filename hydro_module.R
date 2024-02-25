trans_height <- function(){
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
}

layer_water_depth <- function(){
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
}

wat_k_mean <- function(){
  for (x in 1:x_extent){
    for (y in 1:y_extent){
      if (activation_status[x,y]=="on"){
        if (water_table[x,y] >= transmissivity[x,y,no_layers-1,1]){
          wk_mean[x,y] <- transmissivity[x,y,no_layers-1,2]
          if (wk_mean[x,y] > transmax) {transmax <- wk_mean[x,y]}
          wk_mean[x,y] <- wk_mean[x,y]/transmissivity[x,y,no_layers-1,1]
        }
        else if (water_table[x,y] < transmissivity[x,y,1,1]){
          wk_mean[x,y] <- layer_attributes[x,y,1,2] #
          if (transmissivity[x,y,1,2] > transmax) {transmax <- transmissivity[x,y,1,2]} #
          #
        }
        else {
          for (z in 1:(no_layers-2)){
            if (water_table[x,y] == transmissivity[x,y,z,1]){
              wk_mean[x,y] <- transmissivity[x,y,z,2]
              if (wk_mean[x,y] > transmax) {transmax <- wk_mean[x,y]}
              wk_mean[x,y] <- wk_mean[x,y]/transmissivity[x,y,z,1] #
            }
            else if (water_table[x,y] > transmissivity[x,y,z,1]){
              if (water_table[x,y] < transmissivity[x,y,z+1,1]){
                wk_mean[x,y] <- transmissivity[x,y,z,2]+layer_attributes[x,y,z+1,2]*(water_table[x,y]-transmissivity[x,y,z,1])
                if (wk_mean[x,y] > transmax) {transmax <- wk_mean[x,y]}
                wk_mean[x,y] <- wk_mean[x,y]/water_table[x,y] #
              }
            }
          }
        }
      }
    }
  }
}

move_water <- function(){
  x_flux <- y_flux <- array(NA,dim=c(x_extent,y_extent))

  for (x in 1:(x_extent-1)){
    for (y in 2:(y_extent-1)){
      if ((activation_status[x,y]=="neu")&(activation_status[x+1,y]=="on")){x_flux[x,y] <- 0}
      else if ((activation_status[x,y]=="diri")&(activation_status[x+1,y]=="on")){
        x_flux[x,y] <- wk_mean[x+1,y]*((water_table[x,y]+water_table[x+1,y])/2)*
          (base_altitude[x,y]+water_table[x,y]-base_altitude[x+1,y]-water_table[x+1,y])*2*timestep
      }
      else if ((activation_status[x,y]=="on")&(activation_status[x+1,y]=="neu")){x_flux[x,y] <- 0}
      else if ((activation_status[x,y]=="on")&(activation_status[x+1,y]=="diri")){
        x_flux[x,y] <- wk_mean[x,y]*((water_table[x,y]+water_table[x+1,y])/2)*
          (base_altitude[x,y]+water_table[x,y]-base_altitude[x+1,y]-water_table[x+1,y])*2*timestep
      }
      else if ((activation_status[x,y]=="on")&(activation_status[x+1,y]=="on")){
        x_flux[x,y] <- (2*wk_mean[x,y]*wk_mean[x+1,y])/(wk_mean[x,y]+wk_mean[x+1,y])*((water_table[x,y]+water_table[x+1,y])/2)*
          (base_altitude[x,y]+water_table[x,y]-base_altitude[x+1,y]-water_table[x+1,y])*timestep
      }
    }
  }
  for (y in 1:(y_extent-1)){
    for (x in 2:(x_extent-1)){
      if ((activation_status[x,y]=="neu")&(activation_status[x,y+1]=="on")){y_flux[x,y] <- 0}
      else if ((activation_status[x,y]=="diri")&(activation_status[x,y+1]=="on")){
        y_flux[x,y] <- wk_mean[x,y+1]*((water_table[x,y]+water_table[x,y+1])/2)*
          (base_altitude[x,y]+water_table[x,y]-base_altitude[x,y+1]-water_table[x,y+1])*2*timestep
      }
      else if ((activation_status[x,y]=="on")&(activation_status[x,y+1]=="neu")){y_flux[x,y] <- 0}
      else if ((activation_status[x,y]=="on")&(activation_status[x,y+1]=="diri")){
        y_flux[x,y] <- wk_mean[x,y]*((water_table[x,y]+water_table[x,y+1])/2)*
          (base_altitude[x,y]+water_table[x,y]-base_altitude[x,y+1]-water_table[x,y+1])*2*timestep
      }
      else if ((activation_status[x,y]=="on")&(activation_status[x+1,y]=="on")){
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
}

water_table_update <- function(){
  for (x in 2:(x_extent-1)){
    for (y in 2:(y_extent-1)){
      if (activation_status[x,y]=="on"){
        if (water_table[x,y] <= layer_storage[x,y,1,1]){
          marker <- 1
        }
        else {
          for (z in 1:(no_layers-1)){
            if (water_table[x,y] > layer_storage[x,y,z,1]){
              if (water_table[x,y] <= layer_storage[x,y,z+1,1]){
                marker <- z+1
              }
              else if (water_table[x,y] > layer_storage[x,y,no_layers,1]){
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
                water_table[x,y] <- water_table[x,y]+water_change[x,y]/layer_attributes(x,y,z,3)
                break
              }
              else {
                if (z==no_layers){
                  water_table[x,y] <- layer_storage[x,y,z,1]
                  break
                }
                else {
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
              }
              else {
                test_depth <- layer_storage[x,y,z,2]*(water_table[x,y]-layer_storage[x,y,z-1,1])/layer_attributes[x,y,z,1]
              }
              if (abs(water_change[x,y]) <= test_depth){
                water_table[x,y] <- water_table[x,y]+water_change[x,y]/layer_attributes[x,y,z,3]
                break
              }
              else {
                if (z==1){
                  water_table[x,y] <- 0
                  break
                }
                else {
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
}