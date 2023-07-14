#include "Uwb_Location/trilateration.h"
#include "Uwb_Location/uwb.h"
#include "std_msgs/String.h" //ros定义的String数据类型
#include <iostream>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <serial/serial.h>
#include <string.h>
#include <string>

using namespace std;
unsigned char receive_buf[200] = { 0 };
vec3d report;
Quaternion q;
int result = 0;
float velocityac[3], angleac[3];
Quaternion Q;

#define baseHeight 2.10
#define MAX_DATA_NUM 1024 // 传消息内容最大长度
#define DataHead 'm'
#define DataTail '\n'
unsigned char BufDataFromCtrl[MAX_DATA_NUM];
int BufCtrlPosit_w = 0, BufCtrlPosit_r = 0;
int DataRecord = 0, rcvsign = 0;
int range[8] = { -1 };
vec3d anchorArray[8];

void receive_deal_func()
{
    for (int i = 0; i < 8; i++)
    {
        range[i] = -1;
    }
    if ((receive_buf[0] == 'm') && (receive_buf[1] == 'c')) {
        int aid, tid, lnum, seq, mask;
        int rangetime;
        char role;
        int data_len = strlen((char*)receive_buf);
        // printf("lenmc = %d\n", data_len);
        if (data_len == 106) {
            int n = sscanf((char*)receive_buf, "mc %x %x %x %x %x %x %x %x %x %x %x %x %c%d:%d",
                &mask,
                &range[0], &range[1], &range[2], &range[3],
                &range[4], &range[5], &range[6], &range[7],
                &lnum, &seq, &rangetime, &role, &tid, &aid);
            printf("mask=0x%02x\nrange[0]=%d(mm)\nrange[1]=%d(mm)\nrange[2]=%d(mm)\nrange[3]=%d(mm)\nrange[4]=%d(mm)\nrange[5]=%d(mm)\nrange[6]=%d(mm)\nrange[7]=%d(mm)\r\n",
                mask, range[0], range[1], range[2], range[3],
                range[4], range[5], range[6], range[7]);
        } else if (data_len == 70) {
            int n = sscanf((char*)receive_buf, "mc %x %x %x %x %x %x %x %x %c%d:%d",
                &mask,
                &range[0], &range[1], &range[2], &range[3],
                &lnum, &seq, &rangetime, &role, &tid, &aid);
            printf("mask=0x%02x\nrange[0]=%d(mm)\nrange[1]=%d(mm)\nrange[2]=%d(mm)\nrange[3]=%d(mm)\r\n",
                mask, range[0], range[1], range[2], range[3]);
        } else {
            return;
        }
    } else if ((receive_buf[0] == 'm') && (receive_buf[1] == 'i')) {
        float rangetime2;
        float acc[3], gyro[3];

        // mi,981.937,0.63,NULL,NULL,NULL,-2.777783,1.655664,9.075048,-0.004788,-0.014364,-0.001596,T0     //13
        // mi,3.710,0.55,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1.327881,0.653174,9.577490,-0.004788,-0.013300,-0.002128,T0    //17
        char *ptr, *retptr;
        ptr = (char*)receive_buf;
        char cut_data[30][12];
        int cut_count = 0;

        while ((retptr = strtok(ptr, ",")) != NULL) {
            // printf("%s\n", retptr);
            strcpy(cut_data[cut_count], retptr);
            ptr = NULL;
            cut_count++;
            if (cut_count >= 29)
                break;
        }

        rangetime2 = atof(cut_data[1]);

        if (cut_count == 13) // 4anchors
        {
            for (int i = 0; i < 4; i++) {
                if (strcmp(cut_data[i + 2], "NULL")) {
                    range[i] = atof(cut_data[i + 2]) * 1000;
                } else {
                    range[i] = -1;
                }
            }

            for (int i = 0; i < 3; i++) {
                acc[i] = atof(cut_data[i + 6]);
            }

            for (int i = 0; i < 3; i++) {
                gyro[i] = atof(cut_data[i + 9]);
            }

            printf("rangetime = %.3f\n", rangetime2);
            printf("range[0] = %d\n", range[0]);
            printf("range[1] = %d\n", range[1]);
            printf("range[2] = %d\n", range[2]);
            printf("range[3] = %d\n", range[3]);
            printf("acc[0] = %.3f\n", acc[0]);
            printf("acc[1] = %.3f\n", acc[1]);
            printf("acc[2] = %.3f\n", acc[2]);
            printf("gyro[0] = %.3f\n", gyro[0]);
            printf("gyro[1] = %.3f\n", gyro[1]);
            printf("gyro[2] = %.3f\n", gyro[2]);
        } else if (cut_count == 17) // 8anchors
        {
            for (int i = 0; i < 8; i++) {
                if (strcmp(cut_data[i + 2], "NULL")) {
                    range[i] = atof(cut_data[i + 2]) * 1000;
                } else {
                    range[i] = -1;
                }
            }

            for (int i = 0; i < 3; i++) {
                acc[i] = atof(cut_data[i + 6 + 4]);
            }

            for (int i = 0; i < 3; i++) {
                gyro[i] = atof(cut_data[i + 9 + 4]);
            }

            printf("rangetime = %.3f\n", rangetime2);
            printf("range[0] = %d\n", range[0]);
            printf("range[1] = %d\n", range[1]);
            printf("range[2] = %d\n", range[2]);
            printf("range[3] = %d\n", range[3]);
            printf("range[4] = %d\n", range[4]);
            printf("range[5] = %d\n", range[5]);
            printf("range[6] = %d\n", range[6]);
            printf("range[7] = %d\n", range[7]);
            printf("acc[0] = %.3f\n", acc[0]);
            printf("acc[1] = %.3f\n", acc[1]);
            printf("acc[2] = %.3f\n", acc[2]);
            printf("gyro[0] = %.3f\n", gyro[0]);
            printf("gyro[1] = %.3f\n", gyro[1]);
            printf("gyro[2] = %.3f\n", gyro[2]);
        } else {
            return;
        }
    } else {
        puts("no range message");
        return;
    }

    // A0 uint:m
    anchorArray[0].x = 0.0;
    anchorArray[0].y = 0.0;
    anchorArray[0].z = baseHeight;
    // A1 uint:m
    anchorArray[1].x = 5.0;
    anchorArray[1].y = 0.0;
    anchorArray[1].z = baseHeight;
    // A2 uint:m
    anchorArray[2].x = 0.0;
    anchorArray[2].y = 5.0;
    anchorArray[2].z = baseHeight;
    // A3 uint:m
    anchorArray[3].x = 5.0;
    anchorArray[3].y = 5.0;
    anchorArray[3].z = baseHeight;
    // // A4 uint:m
    // anchorArray[4].x = 0;
    // anchorArray[4].y = 0;
    // anchorArray[4].z = 0;
    // // A5 uint:m
    // anchorArray[5].x = 0;
    // anchorArray[5].y = 0;
    // anchorArray[5].z = 0;
    // // A6 uint:m
    // anchorArray[6].x = 0;
    // anchorArray[6].y = 0;
    // anchorArray[6].z = 0;
    // // A7 uint:m
    // anchorArray[7].x = 0;
    // anchorArray[7].y = 0;
    // anchorArray[7].z = 0;

    result = GetLocation(&report, &anchorArray[0], &range[0]);

    printf("result = %d\n", result);
    printf("x = %f\n",report.x);
    printf("y = %f\n",report.y);
    printf("z = %f\n",report.z);

}

void CtrlSerDataDeal()
{
    unsigned char middata = 0;
    static unsigned char dataTmp[MAX_DATA_NUM] = { 0 };

    while (BufCtrlPosit_r != BufCtrlPosit_w) {
        middata = BufDataFromCtrl[BufCtrlPosit_r];
        BufCtrlPosit_r = (BufCtrlPosit_r == MAX_DATA_NUM - 1) ? 0 : (BufCtrlPosit_r + 1);

        if (((middata == DataHead)) && (rcvsign == 0)) // 收到头
        {
            rcvsign = 1; // 开始了一个数据帧
            dataTmp[DataRecord++] = middata; // 数据帧接收中
        } else if ((middata != DataTail) && (rcvsign == 1)) {
            dataTmp[DataRecord++] = middata; // 数据帧接收中
        } else if ((middata == DataTail) && (rcvsign == 1)) // 收到尾
        {
            if (DataRecord != 1) {
                rcvsign = 0;
                dataTmp[DataRecord++] = middata;
                dataTmp[DataRecord] = '\0';

                strncpy((char*)receive_buf, (char*)dataTmp, DataRecord);
                printf("receive_buf = %slen = %d\n", receive_buf, DataRecord);
                receive_deal_func(); /*调用处理函数*/
                bzero(receive_buf, sizeof(receive_buf));

                DataRecord = 0;
            }
        }
    }
}

int main(int argc, char** argv)
{
    setlocale(LC_ALL, "");
    std_msgs::String msg;
    std_msgs::String msg_mc;
    int data_size;
    int n;
    int cnt = 0;
    ros::init(argc, argv, "uwb_imu_node"); // 发布imu,uwb节点
    // 创建句柄（虽然后面没用到这个句柄，但如果不创建，运行时进程会出错）
    ros::NodeHandle nh;
    ros::NodeHandle nh1;
    ros::Publisher uwb_publisher = nh.advertise<Uwb_Location::uwb>("/uwb/data", 1000); // 发布uwb数据  话题名 队列大小
    ros::Publisher IMU_read_pub = nh.advertise<sensor_msgs::Imu>("imu/data", 1000); // 发布imu话题

    // 创建一个serial类
    serial::Serial sp;
    // 创建timeout
    serial::Timeout to = serial::Timeout::simpleTimeout(11);
    // 设置要打开的串口名称
    sp.setPort("/dev/ttyUSB0");
    // 设置串口通信的波特率
    sp.setBaudrate(115200);
    // 串口设置timeout
    sp.setTimeout(to);

    try {
        // 打开串口
        sp.open();
    } catch (serial::IOException& e) {
        ROS_ERROR_STREAM("Unable to open port.");
        return -1;
    }

    // 判断串口是否打开成功
    if (sp.isOpen()) {
        ROS_INFO_STREAM("/dev/ttyUSB0 is opened.");
    } else {
        return -1;
    }

    // ros::Rate loop_rate(11);

    // 发布uwb话题
    Uwb_Location::uwb uwb_data;
    // 打包IMU数据
    sensor_msgs::Imu imu_data;

    while (ros::ok()) {
        // 获取缓冲区内的字节数
        size_t len = sp.available();
        if (len > 0) {
            unsigned char usart_buf[1024] = { 0 };
            sp.read(usart_buf, len);
            // printf("uart_data = %s\n", usart_buf);

            unsigned char* pbuf;
            unsigned char buf[2014] = { 0 };

            pbuf = (unsigned char*)usart_buf;
            memcpy(&buf[0], pbuf, len);

            int reallength = len;
            int i;
            if (reallength != 0) {
                for (i = 0; i < reallength; i++) {
                    BufDataFromCtrl[BufCtrlPosit_w] = buf[i];

                    BufCtrlPosit_w = (BufCtrlPosit_w == (MAX_DATA_NUM - 1)) ? 0 : (1 + BufCtrlPosit_w);
                }
            }
            CtrlSerDataDeal();
            //---------------------------------UWB----------------------------------------------------
            uwb_data.time = ros::Time::now();
            uwb_data.x = report.x;
            uwb_data.y = report.y;

            // // 纠正z的算法
            // double proj[4] = {
            //     sqrt(report.x * report.x + report.y * report.y),
            //     sqrt((report.x - anchorArray[1].x) * (report.x - anchorArray[1].x) + (report.y - anchorArray[1].y) * (report.y - anchorArray[1].y)),
            //     sqrt((report.x - anchorArray[2].x) * (report.x - anchorArray[2].x) + (report.y - anchorArray[2].y) * (report.y - anchorArray[2].y)),
            //     sqrt((report.x - anchorArray[3].x) * (report.x - anchorArray[3].x) + (report.y - anchorArray[3].y) * (report.y - anchorArray[3].y))
            // };

            // double h0 = sqrt((double)(range[0] * range[0]) / 1000000.00 - proj[0] * proj[0] > 0 ? (double)(range[0] * range[0]) / 1000000.00 - proj[0] * proj[0] : 0);
            // double h1 = sqrt((double)(range[1] * range[1]) / 1000000.00 - proj[1] * proj[1] > 0 ? (double)(range[1] * range[1]) / 1000000.00 - proj[1] * proj[1] : 0);
            // double h2 = sqrt((double)(range[2] * range[2]) / 1000000.00 - proj[2] * proj[2] > 0 ? (double)(range[2] * range[2]) / 1000000.00 - proj[2] * proj[2] : 0);

            // double h;
            // if (range[3] != 0) {
            //     double h3 = sqrt((double)(range[3] * range[3]) / 1000000.00 - proj[3] * proj[3] > 0 ? (double)(range[3] * range[3]) / 1000000.00 - proj[3] * proj[3] : 0);
            //     // printf("h0 = %f, h1 = %f, h2 = %f, h3 = %f\n", h0, h1, h2, h3);
            //     h = (h0 + h1 + h2 + h3) / 4;
            // } else
            //     h = (h0 + h1 + h2) / 3;

            // uwb_data.z = h + baseHeight;

            // printf("x = %f\n", uwb_data.x);
            // printf("y = %f\n", uwb_data.y);
            // printf("z = %f\n", uwb_data.z);
            uwb_data.z = report.z;
            // printf("tag.x=%.3f\r\ntag.y=%.3f\r\ntag.z=%.3f\r\n",uwb_data.x,uwb_data.y,uwb_data.z);

            //--------------------------------------IMU------------------------------------------------
            // 发布imu话题
            imu_data.header.stamp = uwb_data.time;
            imu_data.header.frame_id = "base_link";
            imu_data.linear_acceleration.x = velocityac[0];
            imu_data.linear_acceleration.y = velocityac[1];
            imu_data.linear_acceleration.z = velocityac[2];

            // 角速度
            imu_data.angular_velocity.x = angleac[0];
            imu_data.angular_velocity.y = angleac[1];
            imu_data.angular_velocity.z = angleac[2];

            // 四元数
            imu_data.orientation.x = Q.q1;
            imu_data.orientation.y = Q.q2;
            imu_data.orientation.z = Q.q3;
            imu_data.orientation.w = Q.q0;

            //--------------------------------------话题发布------------------------------------
            uwb_publisher.publish(uwb_data);
            IMU_read_pub.publish(imu_data);
        }

        ros::spinOnce();
        // loop_rate.sleep();
    }
    // 关闭串口
    sp.close();
    return 0;
}
