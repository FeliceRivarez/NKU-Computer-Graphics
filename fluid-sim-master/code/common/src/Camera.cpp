#include "Camera.h"
#include <iostream>

namespace Glb {

    void Camera::ProcessMove(glm::vec2 offset) {
        mPosition -= offset.x * mSensitiveX * mRight;
        mPosition += offset.y * mSensitiveY * mUp;
        UpdateView();
    }

    void Camera::ProcessRotate(glm::vec2 offset) {
        mYaw = std::fmodf(mYaw - mSensitiveYaw * offset.x, 360.0f);
        mPitch = glm::clamp(mPitch + mSensitivePitch * offset.y, -89.9f, 89.9f);
        UpdateView();
    }

    void Camera::ProcessScale(float offset) {
        mPosition += offset * mSensitiveFront * mFront;
        UpdateView();
    }

    glm::mat4 Camera::GetView() {
        return glm::lookAt(mPosition, mPosition + mFront, mUp);
    }

    glm::mat4 Camera::GetProjection() {
        return glm::perspective(glm::radians(fovyDeg), aspect, nearPlane, farPlane);
    }

    glm::vec3 Camera::GetUp() {
        return mUp;
    }

    glm::vec3 Camera::GetRight() {
        return mRight;
    }

    glm::vec3 Camera::GetFront() {
        return mFront;
    }

    glm::vec3 Camera::GetPosition() {
        return mPosition;
    }

    void Camera::UpdateView() {
        // 更新三个向量
        mFront.x = std::cos(glm::radians(mPitch)) * std::cos(glm::radians(mYaw));
        mFront.y = std::cos(glm::radians(mPitch)) * std::sin(glm::radians(mYaw));
        mFront.z = std::sin(glm::radians(mPitch));
        mFront = -glm::normalize(mFront);

        mRight = glm::normalize(glm::cross(mFront, mWorldUp));
        mUp = glm::normalize(glm::cross(mRight, mFront));
    }
}